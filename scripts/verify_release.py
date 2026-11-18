#!/usr/bin/env python3
"""Verify release contents and exercise each distributable outside the checkout.

The wheel is installed with pip into a temporary target. Zipapps run in isolated
Python mode, and the portable launcher runs after extraction to a new directory.
No scientific engine or network connection is needed for these checks.
"""
from __future__ import annotations

import argparse
import ast
import email.parser
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import sys
import tarfile
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from av_version import VERSION
from checksums import release_artifacts
from release_files import EXCLUDED_PARTS, manual_files, portable_files, runtime_files


class ReleaseError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ReleaseError(message)


def archive_path(name):
    path = PurePosixPath(name)
    require(not path.is_absolute() and '..' not in path.parts and '\\' not in name,
            'Unsafe archive path: ' + name)
    require(not any(part in EXCLUDED_PARTS for part in path.parts),
            'Excluded development resource in archive: ' + name)
    require(path.suffix not in ('.pyc', '.pyo'), 'Bytecode cache in archive: ' + name)
    return path


def source_version(content):
    for node in ast.parse(content.decode('utf-8')).body:
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == 'VERSION'
                                               for target in node.targets):
            return ast.literal_eval(node.value)
    raise ReleaseError('av_version.py has no literal VERSION assignment.')


def check_zip(path, kind):
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)), f'{path.name}: duplicate member names')
        for member in archive.infolist():
            archive_path(member.filename)
            require(not stat.S_ISLNK(member.external_attr >> 16),
                    f'{path.name}: archive contains a link: {member.filename}')
        prefix = 'AlignView/' if kind == 'portable' else ''
        if kind == 'portable':
            require(all(name.startswith(prefix) for name in names), 'Portable ZIP must contain one AlignView/ directory.')
            require(prefix + 'alignview.sh' in names, 'Portable launcher missing.')
            require(prefix + 'alignview.pyz' in names, 'Portable zipapp missing.')
            require(archive.read(prefix + 'alignview.pyz') == (ROOT / 'alignview.pyz').read_bytes(),
                    'Portable archive embeds an outdated zipapp.')
            expected_files = {prefix + file.relative_to(ROOT).as_posix() for file in portable_files(ROOT)}
            expected_files.add(prefix + 'alignview.pyz')
            require(set(names) == expected_files, 'Portable ZIP differs from its declared file inventory.')
            for file in portable_files(ROOT):
                name = prefix + file.relative_to(ROOT).as_posix()
                require(archive.read(name) == file.read_bytes(), f'{path.name}: stale resource {name}')
        require(source_version(archive.read(prefix + 'av_version.py')) == VERSION,
                f'{path.name}: embedded source version differs from {VERSION}')
        for module in runtime_files(ROOT):
            name = prefix + module.name
            require(name in names, f'{path.name}: missing {name}')
            require(archive.read(name) == module.read_bytes(), f'{path.name}: stale runtime module {module.name}')
        runtime_members = {name for name in names if name.startswith(prefix)
                           and '/' not in name[len(prefix):] and name.endswith('.py')}
        expected = {prefix + module.name for module in runtime_files(ROOT)}
        if kind == 'zipapp':
            expected.add('__main__.py')
        require(runtime_members == expected, f'{path.name}: undeclared root Python modules')
        manual_prefix = 'docs/' if kind == 'portable' else 'alignview_docs/'
        expected_manual = set()
        for document in manual_files(ROOT):
            name = prefix + manual_prefix + document.relative_to(ROOT / 'docs').as_posix()
            expected_manual.add(name)
            require(name in names, f'{path.name}: missing offline guide {document.name}')
            require(archive.read(name) == document.read_bytes(), f'{path.name}: stale offline resource {document.name}')
        if kind == 'zipapp':
            require(set(names) == expected | expected_manual, 'Zipapp contains undeclared resources.')
        if kind == 'wheel':
            metadata_name = next((name for name in names if name.endswith('.dist-info/METADATA')), '')
            require(bool(metadata_name), 'Wheel metadata missing.')
            metadata = email.parser.Parser().parsestr(archive.read(metadata_name).decode('utf-8'))
            require(metadata['Name'] == 'alignview' and metadata['Version'] == VERSION,
                    'Wheel metadata version/name mismatch.')
            entrypoints = metadata_name.rsplit('/', 1)[0] + '/entry_points.txt'
            require('alignview = alignview:main' in archive.read(entrypoints).decode('utf-8'),
                    'Installed alignview console entry point missing.')
            metadata_prefix = metadata_name.rsplit('/', 1)[0] + '/'
            require(all(name in expected | expected_manual or name.startswith(metadata_prefix) for name in names),
                    'Wheel contains undeclared application resources.')


def check_source(path):
    with tarfile.open(path, 'r:gz') as archive:
        members = archive.getmembers()
        prefix = f'alignview-{VERSION}/'
        names = {member.name for member in members}
        require(len(members) == len(names), 'Source archive contains duplicate members.')
        for member in members:
            archive_path(member.name)
            require(member.name == prefix[:-1] or member.name.startswith(prefix),
                    'Source archive has an unexpected top-level directory.')
            require(not member.issym() and not member.islnk(), 'Source archive contains a link: ' + member.name)
            require(member.isfile() or member.isdir(), 'Source archive contains a special file: ' + member.name)
        for module in runtime_files(ROOT):
            name = prefix + module.name
            require(name in names, 'Source archive missing ' + module.name)
            stream = archive.extractfile(name)
            require(stream is not None and stream.read() == module.read_bytes(),
                    'Source archive contains stale module ' + module.name)
        runtime_names = {name for name in names if name.startswith(prefix)
                         and '/' not in name[len(prefix):] and name.endswith('.py')}
        require(runtime_names == {prefix + module.name for module in runtime_files(ROOT)},
                'Source archive contains undeclared root Python modules.')
        for document in manual_files(ROOT):
            name = prefix + 'docs/' + document.relative_to(ROOT / 'docs').as_posix()
            require(name in names, 'Source archive missing offline resource ' + document.name)
            stream = archive.extractfile(name)
            require(stream is not None and stream.read() == document.read_bytes(),
                    'Source archive contains stale offline resource ' + document.name)
        metadata_file = archive.extractfile(prefix + 'PKG-INFO')
        require(metadata_file is not None, 'Source archive metadata missing.')
        metadata = email.parser.Parser().parsestr(metadata_file.read().decode('utf-8'))
        require(metadata['Name'] == 'alignview' and metadata['Version'] == VERSION,
                'Source archive metadata version/name mismatch.')


def check_checksums(paths):
    manifest = ROOT / 'dist/SHA256SUMS'
    require(manifest.is_file(), 'Run scripts/checksums.py first.')
    entries = {}
    for line in manifest.read_text(encoding='utf-8').splitlines():
        digest, name = line.split('  ', 1)
        require(Path(name).name == name and name not in entries, 'Invalid checksum entry: ' + name)
        require(re.fullmatch('[0-9a-f]{64}', digest), 'Invalid SHA-256 digest for ' + name)
        entries[name] = digest
    require(set(entries) == {path.name for path in paths}, 'Checksums must cover exactly the current release set.')
    for path in paths:
        require(hashlib.sha256(path.read_bytes()).hexdigest() == entries[path.name],
                'Checksum mismatch: ' + path.name)


def run(command, cwd, *, env=None):
    result = subprocess.run(command, cwd=cwd, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, timeout=90)
    require(result.returncode == 0, 'Command failed: ' + ' '.join(map(str, command)) + '\n' + result.stdout[-6000:])
    return result.stdout


def smoke(command, workspace, label):
    caller = workspace / label
    caller.mkdir()
    env = dict(os.environ)
    env['NO_COLOR'] = '1'
    env['ALIGNVIEW_LANG'] = 'fr'  # Old environment settings cannot change the English interface.
    version = run([*command, '--version'], caller, env=env)
    require(VERSION in version, label + ': version output differs from release version')
    help_text = run([*command, '--help'], caller, env=env)
    require('--help' in help_text and '--lang' not in help_text, label + ': invalid English help output')
    manual = run([*command, '--docs', 'quickstart'], caller, env=env)
    require(manual.startswith('# Getting started'), label + ': bundled manual unavailable')
    run([*command, '--demo', '--no-menu'], caller, env=env)
    (caller / 'input.fasta').write_text('>alpha\nACGT-AC\n>beta\nACGTAAC\n>gamma\nACGTTAC\n', encoding='utf-8')
    run([*command, 'view', 'input.fasta', '--no-menu'], caller, env=env)
    run([*command, 'report', 'input.fasta', '--html-out', 'reports/result.html'], caller, env=env)
    report = caller / 'reports/result.html'
    require(report.is_file(), label + ': caller-relative HTML output missing')
    markup = report.read_text(encoding='utf-8')
    require('lang="en"' in markup or "lang='en'" in markup, label + ': report language is not English')
    require('alpha' in markup and 'gamma' in markup, label + ': report lost sequence data')
    require(not (caller / 'session.json').exists(), label + ': view/export unexpectedly saved a project')
    run([*command, 'input.fasta', '--msa', 'view', '--save', '--out', 'saved project', '--no-menu'], caller, env=env)
    require((caller / 'saved project/session.json').is_file(), label + ': caller-relative project missing')
    run([*command, 'resume', 'saved project', '--no-menu'], caller, env=env)
    run([*command, 'report', 'saved project', '--html-out', 'reports/resumed.html'], caller, env=env)
    require((caller / 'reports/resumed.html').is_file(), label + ': project report command failed')
    print(label + ': version, help, manual, demo, reports and project resume passed')


def verify_release(*, contents_only=False):
    paths = release_artifacts(ROOT)
    check_checksums(paths)
    by_suffix = {('wheel' if path.suffix == '.whl' else 'portable' if path.suffix == '.zip' else
                  'zipapp' if path.suffix == '.pyz' else 'source'): path for path in paths}
    for kind in ('wheel', 'portable', 'zipapp'):
        check_zip(by_suffix[kind], kind)
    check_source(by_suffix['source'])
    require((ROOT / 'alignview.pyz').read_bytes() == by_suffix['zipapp'].read_bytes(),
            'Developer zipapp differs from the versioned release zipapp.')
    print(f'AlignView {VERSION}: release inventory, source freshness, metadata and checksums passed')
    if contents_only:
        return
    with tempfile.TemporaryDirectory(prefix='alignview_release_check_') as temp:
        workspace = Path(temp)
        target = workspace / 'installed'
        # pip can come from the build environment; application imports below are isolated.
        run([sys.executable, '-m', 'pip', 'install', '--no-index', '--no-deps', '--target', str(target),
             str(by_suffix['wheel'])], workspace)
        bootstrap = ('import sys;sys.path.insert(0,' + repr(str(target)) + ');'
                     'from alignview import main;raise SystemExit(main())')
        smoke([sys.executable, '-I', '-c', bootstrap], workspace, 'installed-wheel')
        smoke([sys.executable, '-I', str(by_suffix['zipapp'])], workspace, 'isolated-zipapp')
        with zipfile.ZipFile(by_suffix['portable']) as archive:
            archive.extractall(workspace / 'portable')
        launcher = workspace / 'portable/AlignView/alignview.sh'
        smoke(['bash', str(launcher)], workspace, 'portable-launcher')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contents-only', action='store_true', help='Check all payloads/checksums without running programs')
    try:
        verify_release(contents_only=parser.parse_args().contents_only)
    except (ValueError, OSError, subprocess.TimeoutExpired, zipfile.BadZipFile, tarfile.TarError) as exc:
        raise SystemExit(str(exc)) from exc
