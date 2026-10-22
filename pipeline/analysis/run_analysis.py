"""Advanced JSON function interface; marked paths are relative to the job file."""
from __future__ import annotations
import argparse
import importlib
import json
from pathlib import Path

MODULES = {'fig2_insample', 'fig2_spatial', 'gene_characterization', 'gene_groups',
           'fig3_hvg', 'fig3_heatmaps', 'fig3_summary', 'fig4_fov',
           'fig4_fov_boxplots', 'fig4_adjacent', 'fig5_resolution', 'fig6_usability', 'gene_selection'}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job', required=True, type=Path)
    parser.add_argument('--result-csv', type=Path)
    args = parser.parse_args(argv)
    job = json.loads(args.job.read_text(encoding='utf-8'))
    if job['module'] not in MODULES:
        parser.error(f"Unknown module: {job['module']}")
    base = args.job.resolve().parent
    def decode(value):
        if isinstance(value, dict):
            if set(value) in ({'__path__'}, {'__csv__'}, {'__image__'}, {'__npy__'}):
                key, val = next(iter(value.items()))
                path = str((base / val).resolve())
                if key == '__path__':
                    return path
                if key == '__csv__':
                    import pandas as pd
                    return pd.read_csv(path)
                import numpy as np
                if key == '__npy__':
                    return np.load(path, allow_pickle=False)
                from PIL import Image
                return np.asarray(Image.open(path))
            return {k: decode(v) for k, v in value.items()}
        if isinstance(value, list):
            return [decode(v) for v in value]
        return value
    module = importlib.import_module('pipeline.analysis.' + job['module'])
    result = getattr(module, job['function'])(**decode(job.get('kwargs', {})))
    if args.result_csv:
        import pandas as pd
        if not isinstance(result, pd.DataFrame):
            raise TypeError('--result-csv requires a DataFrame return value')
        args.result_csv.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.result_csv, index=False)
    print(f"Completed {job['module']}.{job['function']}")


if __name__ == '__main__':
    main()
