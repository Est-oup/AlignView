#!/usr/bin/env python3
"""Prepare a checked Bioconda recipe from an exact AlignView source archive.

Local mode produces a buildable candidate, never a submission claim. Public mode
fetches an HTTPS archive and verifies its SHA256 and runtime/manual contents
against the current checkout before writing a submission directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path,PurePosixPath
import re
import shutil
import sys
import tarfile
import tempfile
import urllib.parse
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_version import VERSION
from release_files import runtime_files,manual_files

MAX_ARCHIVE_BYTES=64*1024*1024


def source_inventory(root=ROOT):
    root=Path(root)
    return [*runtime_files(root),*manual_files(root),root/'pyproject.toml',root/'LICENSE',root/'README.md']


def verify_archive(archive,root=ROOT):
    """Inspect without extracting or executing downloaded upstream code."""
    root=Path(root)
    if Path(archive).stat().st_size>MAX_ARCHIVE_BYTES:
        raise ValueError('The source archive exceeds the 64 MiB preparation limit.')
    with tarfile.open(archive,'r:*') as source:
        members={}
        prefixes=set()
        for entry in source.getmembers():
            path=PurePosixPath(entry.name)
            if path.is_absolute() or '..' in path.parts or '\\' in entry.name:
                raise ValueError('Unsafe source archive member: '+entry.name)
            if not path.parts:continue
            prefixes.add(path.parts[0])
            if entry.issym() or entry.islnk() or not (entry.isfile() or entry.isdir()):
                raise ValueError('Source archives must contain only regular files and directories.')
            if entry.name in members:raise ValueError('Duplicate source archive member: '+entry.name)
            members[entry.name]=entry
        if len(prefixes)!=1:raise ValueError('The source archive must have one top-level directory.')
        prefix=next(iter(prefixes))
        for expected in source_inventory(root):
            name=prefix+'/'+expected.relative_to(root).as_posix()
            entry=members.get(name)
            if entry is None or not entry.isfile():raise ValueError('Source archive is missing '+expected.relative_to(root).as_posix())
            stream=source.extractfile(entry)
            if stream is None or entry.size!=expected.stat().st_size or stream.read()!=expected.read_bytes():
                raise ValueError('Source archive differs from current sources: '+expected.relative_to(root).as_posix())
    return hashlib.sha256(Path(archive).read_bytes()).hexdigest()


def validate_public_url(url,version=VERSION):
    parts=urllib.parse.urlsplit(url)
    if parts.scheme!='https' or not parts.hostname or parts.username or parts.password or parts.fragment or parts.query:
        raise ValueError('Use a public HTTPS source archive URL without credentials, a query string or a fragment.')
    if any(c.isspace() for c in url):raise ValueError('Source URLs cannot contain whitespace.')
    if not parts.path.endswith(('.tar.gz','.tgz','.tar.bz2','.tar.xz')):
        raise ValueError('Use a stable versioned source tarball URL, not a branch or repository URL.')
    if version not in urllib.parse.unquote(parts.path):
        raise ValueError('The source URL must identify the current release version '+version+'.')
    if re.search(r'/(?:heads|branches)/',parts.path):
        raise ValueError('A moving branch is not a release source URL.')
    return url


def download(url,target):
    request=urllib.request.Request(url,headers={'User-Agent':'AlignView-Bioconda-preparation/'+VERSION})
    with urllib.request.urlopen(request,timeout=45) as response,target.open('wb') as output:
        if urllib.parse.urlsplit(response.url).scheme!='https':raise ValueError('The source download redirected away from HTTPS.')
        total=0
        while block:=response.read(1024*1024):
            total+=len(block)
            if total>MAX_ARCHIVE_BYTES:raise ValueError('The source download exceeds 64 MiB.')
            output.write(block)


def render_recipe(url,sha256,maintainer,root=ROOT):
    if not re.fullmatch(r'[0-9a-f]{64}',sha256):raise ValueError('Use the exact 64-character source SHA256.')
    if not re.fullmatch(r'[A-Za-z0-9](?:[A-Za-z0-9-]{0,37}[A-Za-z0-9])?',maintainer):
        raise ValueError('Provide the real GitHub handle of the consenting recipe maintainer.')
    replacements={'@@VERSION@@':VERSION,'@@SOURCE_URL@@':json.dumps(url),
                  '@@SOURCE_SHA256@@':sha256,'@@MAINTAINER@@':maintainer}
    result=(Path(root)/'packaging/bioconda/meta.yaml.template').read_text(encoding='utf-8')
    for token,value in replacements.items():result=result.replace(token,value)
    if re.search(r'@@[A-Z_]+@@',result):raise ValueError('Unresolved Bioconda recipe placeholder.')
    return result


def prepare(*,archive=None,source_url=None,local=False,expected_sha256=None,maintainer,output,root=ROOT):
    root=Path(root);output=Path(output).expanduser().resolve()
    if output.name!='alignview':raise ValueError('The recipe directory must be named alignview.')
    if output.exists():raise ValueError('The recipe destination already exists; choose a new directory: '+str(output))
    if local:
        if not archive or source_url or expected_sha256:
            raise ValueError('--local requires --archive and cannot use --source-url/--sha256.')
    elif not source_url or not expected_sha256:
        raise ValueError('Submission mode requires --source-url and its independently recorded --sha256.')
    if source_url:validate_public_url(source_url)
    with tempfile.TemporaryDirectory(prefix='alignview_bioconda_prepare_') as temporary:
        if local:
            current=Path(archive).expanduser().resolve()
            if not current.is_file():raise ValueError('Source archive not found: '+str(current))
            url=current.as_uri()
        else:
            current=Path(temporary)/'source.tar.gz'
            download(source_url,current)
            url=source_url
        checksum=verify_archive(current,root)
        if expected_sha256 and checksum!=expected_sha256.lower():raise ValueError('Public source SHA256 does not match --sha256.')
        if archive and not local and hashlib.sha256(Path(archive).read_bytes()).hexdigest()!=checksum:
            raise ValueError('The public source differs from the supplied local archive.')
        recipe=render_recipe(url,checksum,maintainer,root)
        mode='local-candidate' if local else 'public-source-verified'
        evidence=dict(schema=1,package='alignview',version=VERSION,mode=mode,source_url=url,
                      source_sha256=checksum,maintainer=maintainer,
                      verified_files=[p.relative_to(root).as_posix() for p in source_inventory(root)],
                      conda_build_run=False,bioconda_lint_run=False,
                      remaining=(['Publish an immutable versioned source archive, then regenerate without --local.'] if local else [])+
                          ['Run Bioconda lint/build for this exact recipe before submission.','Confirm that the named recipe maintainer agrees to maintain it.'])
        # Stage first so a failed validation never leaves a partial recipe.
        staged=Path(temporary)/'recipe';staged.mkdir()
        (staged/'meta.yaml').write_text(recipe,encoding='utf-8')
        output.parent.mkdir(parents=True,exist_ok=True)
        shutil.copytree(staged,output)
        evidence_path=output.parent/'alignview-preparation.json'
        if evidence_path.exists():
            shutil.rmtree(output)
            raise ValueError('Preparation evidence already exists; use a new parent directory: '+str(evidence_path))
        evidence_path.write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return output,evidence_path


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',type=Path,help='Exact source tarball: mandatory for local mode, optional comparison in public mode')
    parser.add_argument('--source-url',help='Stable versioned HTTPS archive URL for submission mode')
    parser.add_argument('--sha256',help='Independently recorded SHA256 of the public source archive')
    parser.add_argument('--local',action='store_true',help='Prepare a local build candidate; it cannot be submitted to Bioconda')
    parser.add_argument('--maintainer',required=True,help='Real GitHub handle; confirm consent before submission')
    parser.add_argument('--output',type=Path,required=True,help='New recipe directory, ending in recipes/alignview for Bioconda')
    args=parser.parse_args(argv)
    try:
        recipe,evidence=prepare(archive=args.archive,source_url=args.source_url,local=args.local,
            expected_sha256=args.sha256,maintainer=args.maintainer,output=args.output)
    except (ValueError,OSError,tarfile.TarError) as exc:
        parser.exit(1,'Cannot prepare Bioconda recipe: '+str(exc)+'\n')
    print(('Local build candidate' if args.local else 'Public-source recipe')+' prepared: '+str(recipe))
    print('Preparation evidence: '+str(evidence))
    print('Nothing has been published or submitted. Run conda-build and Bioconda lint before submission.')
    return 0


if __name__=='__main__':raise SystemExit(main())
