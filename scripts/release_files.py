"""Shared explicit release inventory and deterministic ZIP metadata.

pyproject.toml declares the runtime modules. Resource roots are enumerated so
scientific engines, development caches and caller results never qualify.
"""
from __future__ import annotations
import ast
import datetime
import os
from pathlib import Path
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
RESOURCE_ROOTS = ('docs', 'examples', 'scripts', 'tests', 'packaging', '.github')
ROOT_FILES = ('README.md', 'LICENSE', 'CHANGELOG.md', 'CONTRIBUTING.md',
              'pyproject.toml', 'MANIFEST.in', 'requirements.txt',
              'environment-minimal.yml', 'environment-full.yml', 'alignview.sh')
EXCLUDED_PARTS = {'__pycache__', '.pytest_cache', '.mypy_cache', '.ruff_cache',
                  '.git', 'tools', 'archives', 'build', 'dist', '.venv'}
EXCLUDED_SUFFIXES = {'.pyc', '.pyo', '.zip', '.whl', '.pyz', '.tar', '.tgz', '.gz', '.bz2', '.xz', '.7z'}
DEFAULT_EPOCH = 946684800


def runtime_modules(root=ROOT):
    """Read the declaration without requiring tomllib on Python 3.10."""
    source = (Path(root) / 'pyproject.toml').read_text(encoding='utf-8')
    section = re.search(r'(?ms)^\[tool\.setuptools\]\s*\n(.*?)(?=^\[|\Z)', source)
    found = re.search(r'(?ms)^py-modules\s*=\s*(\[.*?\])', section[1] if section else '')
    if not found:
        raise ValueError('Declare runtime modules under [tool.setuptools].py-modules.')
    modules = ast.literal_eval(found[1])
    if not isinstance(modules, list) or not modules or any(
            not isinstance(name, str) or not re.fullmatch(r'[A-Za-z_]\w*', name)
            for name in modules):
        raise ValueError('Runtime inventory must be a list of simple module names.')
    if len(modules) != len(set(modules)):
        raise ValueError('The runtime inventory contains duplicate modules.')
    return tuple(modules)


def runtime_files(root=ROOT):
    root = Path(root)
    modules = runtime_modules(root)
    paths = [root / (name + '.py') for name in modules]
    for path in paths:
        if path.is_symlink() and not path.resolve().is_relative_to(root.resolve()):
            raise ValueError('Runtime module points outside the application: ' + path.name)
    missing = [path.name for path in paths if not path.is_file()]
    if missing:
        raise ValueError('Declared runtime modules are missing: ' + ', '.join(missing))
    declared = set(modules)
    for path in paths:
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8'), filename=str(path))):
            imports = ([node.module] if isinstance(node, ast.ImportFrom) and node.module else
                       [alias.name for alias in node.names] if isinstance(node, ast.Import) else [])
            for name in imports:
                base = name.split('.', 1)[0]
                if base.startswith('av_') and base not in declared:
                    raise ValueError(f'{path.name} imports undeclared runtime module {base}.')
    return paths


def _eligible(path, root):
    relative = path.relative_to(root)
    if any(part in EXCLUDED_PARTS or part.endswith('.egg-info') for part in relative.parts):
        return False
    if path.suffix in EXCLUDED_SUFFIXES or path.name == '.DS_Store':
        return False
    if path.is_symlink() and not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Release resource points outside the application: ' + str(relative))
    return path.is_file()


def manual_files(root=ROOT):
    """Enumerate packaged offline guides without following external resources."""
    root = Path(root)
    return sorted(path for path in (root / 'docs').rglob('*')
                  if (path.suffix in ('.md', '.svg', '.json') or path.name == '__init__.py')
                  and _eligible(path, root))


def portable_files(root=ROOT):
    root = Path(root)
    paths = [*runtime_files(root), *(root / name for name in ROOT_FILES)]
    for name in RESOURCE_ROOTS:
        paths.extend(path for path in (root / name).rglob('*') if _eligible(path, root))
    missing = [path.name for path in paths if not path.is_file()]
    if missing:
        raise ValueError('Required release files are missing: ' + ', '.join(missing))
    return sorted(set(paths), key=lambda path: path.relative_to(root).as_posix())


def zip_timestamp():
    try:
        epoch = int(os.environ.get('SOURCE_DATE_EPOCH', DEFAULT_EPOCH))
        maximum = int(datetime.datetime(2107, 12, 31, 23, 59, 58, tzinfo=datetime.timezone.utc).timestamp())
        value = datetime.datetime.fromtimestamp(min(maximum, max(315532800, epoch)), datetime.timezone.utc)
    except (ValueError, OverflowError, OSError) as exc:
        raise ValueError('SOURCE_DATE_EPOCH must be a valid Unix timestamp.') from exc
    return (min(value.year, 2107), value.month, value.day, value.hour, value.minute, value.second)


def write_member(archive, name, data, *, executable=False):
    info = zipfile.ZipInfo(name, zip_timestamp())
    info.create_system = 3
    info.external_attr = ((0o100755 if executable else 0o100644) << 16)
    info.compress_type = zipfile.ZIP_DEFLATED
    archive.writestr(info, data, compresslevel=9)
