"""Execute a copy of a platform notebook with explicit path/setup parameters."""
from __future__ import annotations
import argparse
import ast
import json
from pathlib import Path
import re


def override_python(source, parameters):
    try:
        tree = ast.parse(source)
    except SyntaxError:
        # Preserve IPython shell cells. Their preceding configuration cells are parameterized.
        return source, set()
    lines = source.splitlines(keepends=True)
    changes, used = [], set()
    for node in tree.body:
        target = node.targets[0] if isinstance(node, ast.Assign) and len(node.targets) == 1 else node.target if isinstance(node, ast.AnnAssign) else None
        if isinstance(target, ast.Name) and target.id in parameters:
            changes.append((node.lineno - 1, node.end_lineno, f'{target.id} = {parameters[target.id]!r}\n'))
            used.add(target.id)
    for first, last, replacement in reversed(changes):
        lines[first:last] = [replacement]
    return ''.join(lines), used


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--notebook', type=Path, required=True)
    parser.add_argument('--parameters', type=Path, required=True, help='JSON mapping existing setup variables to values')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cwd', type=Path)
    parser.add_argument('--timeout', type=int, default=3600)
    args = parser.parse_args(argv)
    import nbformat
    from nbclient import NotebookClient
    nb = nbformat.read(args.notebook, as_version=4)
    parameters = json.loads(args.parameters.read_text(encoding='utf-8'))
    language = nb.metadata.kernelspec.language.lower()
    used = set()
    for cell in nb.cells:
        if cell.cell_type != 'code':
            continue
        cell.outputs, cell.execution_count = [], None
        if language == 'r':
            for key, value in parameters.items():
                if not isinstance(value, (str, int, float, bool)):
                    raise TypeError('R notebook parameters must be scalar values')
                pattern = rf'(?m)^{re.escape(key)}\s*(?:<-|=)[^\n]*'
                replacement = f'{key} <- ' + ('TRUE' if value is True else 'FALSE' if value is False else json.dumps(value))
                cell.source, count = re.subn(pattern, lambda _: replacement, cell.source)
                if count:
                    used.add(key)
        else:
            cell.source, matched = override_python(cell.source, parameters)
            used.update(matched)
    if set(parameters) - used:
        raise ValueError(f'Unknown setup parameter names: {sorted(set(parameters) - used)}')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    client = NotebookClient(nb, timeout=args.timeout,
        resources={'metadata': {'path': str((args.cwd or args.notebook.parent).resolve())}})
    try:
        client.execute()
    finally:
        # Retain diagnostics in the derived notebook; never overwrite the source notebook.
        nbformat.write(nb, args.output)


if __name__ == '__main__':
    main()
