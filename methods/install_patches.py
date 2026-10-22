"""Install or verify the supplied study modifications in an upstream source checkout."""
import argparse
from pathlib import Path
import shutil

PATCHES = {
    'iSTAR': {'iSTAR/impute.py': 'impute.py'},
    'scstGCN': {'scstGCN/impute.py': 'impute.py',
                'scstGCN/histological_feature_map.py': 'histological_feature_map.py'},
    'XFuse': {'XFuse/visium.py': 'xfuse/convert/visium.py',
              'XFuse/utility.py': 'xfuse/convert/utility.py'},
}

def verify_patches(method, method_dir):
    root = Path(method_dir).expanduser().resolve()
    patch_root = Path(__file__).resolve().parent / 'patches'
    mismatches = []
    for source, destination in PATCHES[method].items():
        installed = root / destination
        if not installed.is_file() or installed.read_bytes() != (patch_root / source).read_bytes():
            mismatches.append(destination)
    if mismatches:
        raise RuntimeError(f"{method} study patches are not installed: {', '.join(mismatches)}. "
                           "Run methods/install_patches.py with --apply first.")

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', choices=PATCHES, required=True)
    parser.add_argument('--method-dir', type=Path, required=True)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--apply', action='store_true', help='Install patches, preserving each original file once')
    action.add_argument('--verify', action='store_true', help='Require byte-identical installed study files')
    args = parser.parse_args(argv)
    root = args.method_dir.expanduser().resolve()
    patch_root = Path(__file__).resolve().parent / 'patches'
    targets = [(patch_root / src, (root / dst).resolve()) for src, dst in PATCHES[args.method].items()]
    for source, destination in targets:
        if not destination.is_relative_to(root) or not destination.is_file():
            raise FileNotFoundError(f"Expected upstream source file: {destination}")
    if args.apply:
        for source, destination in targets:
            if source.read_bytes() == destination.read_bytes():
                continue
            backup = destination.with_name(destination.name + '.srbench-original')
            if not backup.exists():
                shutil.copy2(destination, backup)
            shutil.copyfile(source, destination)
    if args.apply or args.verify:
        verify_patches(args.method, root)
    for source, destination in targets:
        status = 'installed' if source.read_bytes() == destination.read_bytes() else 'not installed'
        print(f"{args.method}: {destination.relative_to(root)}: {status}")
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
