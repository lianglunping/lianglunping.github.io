#!/usr/bin/env python3
"""Read-only release completeness check. Does not upload, build, repair or approve deployment."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys
from typing import Any

SCHEMA = 'rps.release-lock/1'
CORE_REQUIRED = {
    'VERSION', 'README.md', '.gitignore', 'build.py', 'index.html', 'help.html',
    '404.html', 'release.json', 'src/engine.js', 'src/io.js', 'src/review.js',
    'src/gff.js', 'src/audit.js', 'src/audit_ui.js', 'src/deployment.js',
    *{f'src/app_{n:02d}.js' for n in range(1, 6)},
    'tests/tm_golden.json', 'tests/test_engine.cjs', 'tests/test_audit.cjs',
    'tests/test_v020.cjs', 'tests/test_site_browser.py', 'tests/requirements.txt',
    'scripts/preflight_release.py', 'scripts/check_public.py',
}

def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'Duplicate JSON key: {key}')
        result[key] = value
    return result

def safe_path(root: Path, name: str) -> Path:
    if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+(?:/[A-Za-z0-9_.-]+)*', name):
        raise ValueError('Invalid relative release path')
    parts = PurePosixPath(name).parts
    if any(x in {'.', '..', '.git'} for x in parts) or str(PurePosixPath(name)) != name:
        raise ValueError('Unsafe relative release path')
    path = root
    for part in parts:
        path = path / part
        if path.is_symlink():
            raise ValueError(f'Symlink rejected: {name}')
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Path escapes release root')
    return path

def load_lock(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 2_000_000:
        raise ValueError('Missing, symbolic or oversized release lock')
    lock = json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=unique_object)
    if not isinstance(lock, dict) or lock.get('schema') != SCHEMA:
        raise ValueError('Unsupported release-lock schema')
    if not isinstance(lock.get('version'), str) or not lock['version']:
        raise ValueError('Missing release version')
    files = lock.get('files')
    if not isinstance(files, dict) or not (1 <= len(files) <= 2000):
        raise ValueError('Invalid file inventory')
    missing = sorted(CORE_REQUIRED - set(files))
    if missing:
        raise ValueError('Lock omits mandatory files: ' + ', '.join(missing))
    for name, spec in files.items():
        # Validate spelling independently of local file existence.
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+(?:/[A-Za-z0-9_.-]+)*', name):
            raise ValueError('Invalid inventory path')
        if any(x in {'.', '..', '.git'} for x in name.split('/')):
            raise ValueError('Unsafe inventory path')
        if not isinstance(spec, dict) or set(spec) != {'bytes', 'sha256'}:
            raise ValueError('Invalid file descriptor')
        if type(spec['bytes']) is not int or not (0 <= spec['bytes'] <= 64_000_000):
            raise ValueError('Invalid file size')
        if not isinstance(spec['sha256'], str) or not re.fullmatch(r'[0-9a-f]{64}', spec['sha256']):
            raise ValueError('Invalid SHA-256')
    return lock

def inspect(root: Path, lock_path: Path | None = None) -> dict[str, Any]:
    root = Path(root)
    result: dict[str, Any] = {
        'schema': 'rps.preflight-report/1', 'status': 'INVALID_LOCK',
        'contentComplete': False, 'deploymentApproved': False,
        'scientificValidation': 'NOT_ESTABLISHED_BY_THIS_CHECK',
        'verifiedFiles': 0, 'missing': [], 'mismatched': [], 'unsafe': [],
    }
    try:
        if root.is_symlink() or not root.is_dir():
            raise ValueError('Release root must be a real directory')
        lock = load_lock(lock_path or root / 'RELEASE_LOCK.json')
        result['version'] = lock['version']
        result['expectedFiles'] = len(lock['files'])
        result['lockSha256'] = hashlib.sha256((lock_path or root / 'RELEASE_LOCK.json').read_bytes()).hexdigest()
        for name, spec in sorted(lock['files'].items()):
            try:
                path = safe_path(root, name)
            except ValueError:
                result['unsafe'].append(name)
                continue
            if not path.is_file():
                result['missing'].append(name)
                continue
            if path.stat().st_size != spec['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest() != spec['sha256']:
                result['mismatched'].append(name)
                continue
            result['verifiedFiles'] += 1
        version = safe_path(root, 'VERSION')
        if version.is_file() and version.read_text(encoding='utf-8').strip() != lock['version']:
            result['mismatched'].append('VERSION:lock-version')
        if result['unsafe'] or result['mismatched']:
            result['status'] = 'BLOCKED_INCONSISTENT'
        elif result['missing']:
            result['status'] = 'BLOCKED_INCOMPLETE'
        else:
            result['status'] = 'CONTENT_COMPLETE_NOT_DEPLOYMENT_APPROVAL'
            result['contentComplete'] = True
    except (ValueError, OSError, UnicodeError, json.JSONDecodeError) as exc:
        result['error'] = str(exc)
    return result

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--output', type=Path, help='Optional report path; never a source file.')
    args = parser.parse_args()
    root = args.root.resolve()
    if args.output:
        output = args.output.resolve()
        # Reports may be outside root or in qa/, never replace release inputs/lock.
        if output.is_relative_to(root) and not output.is_relative_to(root / 'qa'):
            parser.error('Report output inside release root must be under qa/')
        if args.output.is_symlink():
            parser.error('Report output must not be a symlink')
    report = inspect(args.root, args.manifest)
    text = json.dumps(report, ensure_ascii=False, indent=2) + '\n'
    print(text, end='')
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding='utf-8')
    return 0 if report['contentComplete'] else (3 if report['status'] == 'INVALID_LOCK' else 2)

if __name__ == '__main__':
    sys.exit(main())
