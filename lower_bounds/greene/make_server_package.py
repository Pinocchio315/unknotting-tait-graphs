#!/usr/bin/env python
"""Assemble a self-contained directory for the correction-term sweep on a computing cluster.

The sweep needs the modules of this directory, the form enumerations of `../owens`, and the frozen input
file.  It does not need the KnotInfo database, SnapPy or the rest of the repository; the interpreter on
the server needs only `numpy`, `sympy` and the standard library.

    python make_server_package.py                       # writes greene_server/ and greene_server.tar.gz
    python make_server_package.py --out /tmp/pkg        # somewhere else
    python make_server_package.py --no-archive          # directory only
"""
from __future__ import annotations
import argparse, hashlib, shutil, subprocess, sys, tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
OWENS = HERE.parent / 'owens'

GREENE = ['greene_dinv.py', 'half_integral.py', 'pin_dinv.py', 'greene_ranks.py',
          'greene_sweep.py', 'slurm_greene.sh', 'README.md']
OWENS_MODULES = ['owens_obstruction.py', 'owens_u3.py', 'owens_u4.py', 'linkform.py']
DATA = ['data/greene_targets.json']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, default=HERE / 'greene_server')
    ap.add_argument('--no-archive', action='store_true')
    args = ap.parse_args()

    missing = [name for name in GREENE if not (HERE / name).is_file()] \
        + [name for name in OWENS_MODULES if not (OWENS / name).is_file()] \
        + [name for name in DATA if not (HERE / name).is_file()]
    if missing:
        sys.exit('missing before packaging: ' + ', '.join(missing)
                 + ('\nrun build_greene_targets.py first' if any(d in missing for d in DATA) else ''))

    if args.out.exists():
        shutil.rmtree(args.out)
    (args.out / 'data').mkdir(parents=True)
    for name in GREENE:
        shutil.copy2(HERE / name, args.out / name)
    for name in OWENS_MODULES:                      # flat, so `import owens_obstruction` just works
        shutil.copy2(OWENS / name, args.out / name)
    for name in DATA:
        shutil.copy2(HERE / name, args.out / name)
    (args.out / 'slurm_greene.sh').chmod(0o755)

    manifest = args.out / 'MANIFEST.sha256'
    lines = []
    for path in sorted(p for p in args.out.rglob('*') if p.is_file()):
        lines.append(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(args.out)}')
    manifest.write_text('\n'.join(lines) + '\n')

    total = sum(p.stat().st_size for p in args.out.rglob('*') if p.is_file())
    print(f'{args.out}  ({len(lines) + 1} files, {total // 1024} KiB)')
    for line in lines:
        print('  ' + line.split('  ', 1)[1])

    if not args.no_archive:
        archive = args.out.with_suffix('.tar.gz')
        with tarfile.open(archive, 'w:gz') as tar:
            tar.add(args.out, arcname=args.out.name)
        print(f'\n{archive}  ({archive.stat().st_size // 1024} KiB)')
        print('\ncopy and submit:')
        print(f'  scp {archive.name} SERVER:~/')
        print(f'  ssh SERVER "tar xzf {archive.name}"')
        print(f'  ssh SERVER "cd {args.out.name} && sbatch slurm_greene.sh"')


if __name__ == '__main__':
    main()
