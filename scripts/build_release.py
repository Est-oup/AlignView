#!/usr/bin/env python3
"""Build and verify the current release without publishing it anywhere."""
from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from av_version import VERSION
from release_files import DEFAULT_EPOCH, runtime_files, write_member


def command(arguments, *, env):
    subprocess.run([sys.executable, *arguments], cwd=ROOT, env=env, check=True)


def normalize_wheel(source, target):
    with zipfile.ZipFile(source) as incoming, zipfile.ZipFile(target, 'w') as outgoing:
        for member in sorted(incoming.infolist(), key=lambda entry: entry.filename):
            mode = (member.external_attr >> 16) & 0o777
            write_member(outgoing, member.filename, incoming.read(member), executable=bool(mode & 0o111))


def normalize_source(source, target, epoch):
    # gzip and tar otherwise retain wall-clock times, host owner names and order.
    with tarfile.open(source, 'r:gz') as incoming, target.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=epoch, compresslevel=9) as compressed:
            with tarfile.open(fileobj=compressed, mode='w', format=tarfile.PAX_FORMAT) as outgoing:
                for member in sorted(incoming.getmembers(), key=lambda entry: entry.name):
                    if member.issym() or member.islnk():
                        raise ValueError('Release source archives must not contain links: ' + member.name)
                    member.uid = member.gid = 0
                    member.uname = member.gname = ''
                    member.mtime = epoch
                    member.mode = 0o755 if member.isdir() or member.mode & 0o111 else 0o644
                    member.pax_headers = {}
                    stream = incoming.extractfile(member) if member.isfile() else None
                    outgoing.addfile(member, stream)


def publish_artifact(source, target):
    """Atomically replace an artifact even when build scratch is another mount."""
    with tempfile.TemporaryDirectory(prefix='.alignview-publish-',dir=target.parent) as temporary:
        staged=Path(temporary)/target.name
        shutil.copyfile(source,staged)
        os.replace(staged,target)


def build_release(*, skip_tests=False, contents_only=False):
    runtime_files(ROOT)
    reference = os.environ.get('GITHUB_REF', '')
    if reference.startswith('refs/tags/') and reference.removeprefix('refs/tags/') != 'v' + VERSION:
        raise ValueError(f'Tag {reference} does not match av_version.VERSION={VERSION}.')
    env = dict(os.environ)
    env.setdefault('SOURCE_DATE_EPOCH', str(DEFAULT_EPOCH))
    epoch = int(env['SOURCE_DATE_EPOCH'])
    if epoch < 0:
        raise ValueError('SOURCE_DATE_EPOCH must not be negative.')
    os.environ['SOURCE_DATE_EPOCH'] = str(epoch)
    if not skip_tests:
        command(['-m', 'unittest', 'discover', '-s', 'tests', '-v'], env=env)
    command(['scripts/generate_cli_docs.py', '--check'], env=env)
    command(['scripts/check_docs.py'], env=env)
    with tempfile.TemporaryDirectory(prefix='alignview_release_build_') as temp:
        stage = Path(temp)
        # A clean temporary source copy prevents setuptools' old build/lib from
        # injecting a deleted module into a fresh wheel.
        from release_files import portable_files
        source = stage / 'source'
        for path in portable_files(ROOT):
            target = source / path.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
        built = stage / 'built'
        subprocess.run([sys.executable, '-m', 'build', '--no-isolation', '--outdir', str(built)],
                       cwd=source, env=env, check=True)
        dist = ROOT / 'dist'
        dist.mkdir(exist_ok=True)
        for artifact in sorted(built.iterdir()):
            normalized = stage / artifact.name
            if artifact.suffix == '.whl':
                normalize_wheel(artifact, normalized)
            elif artifact.name.endswith('.tar.gz'):
                normalize_source(artifact, normalized, epoch)
            else:
                raise ValueError('Unexpected build output: ' + artifact.name)
            publish_artifact(normalized, dist / artifact.name)
    command(['scripts/build_pyz.py'], env=env)
    command(['scripts/build_portable.py'], env=env)
    command(['scripts/checksums.py'], env=env)
    command(['-m', 'twine', 'check', str(ROOT / 'dist' / f'alignview-{VERSION}-py3-none-any.whl'),
             str(ROOT / 'dist' / f'alignview-{VERSION}.tar.gz')], env=env)
    command(['scripts/verify_release.py', *(['--contents-only'] if contents_only else [])], env=env)
    print(f'AlignView {VERSION}: release built and verified under {ROOT / "dist"}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-tests', action='store_true', help='Use only after the full suite has passed on these sources')
    parser.add_argument('--contents-only', action='store_true', help='Omit installed/isolated smoke tests when explicitly requested')
    args = parser.parse_args()
    try:
        build_release(skip_tests=args.skip_tests, contents_only=args.contents_only)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit(str(exc)) from exc
