"""Convert trusted, locally generated XFuse tensors to evaluation pickle maps."""
import argparse
import gzip
from pathlib import Path
import pickle


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    import numpy as np
    import torch
    files = sorted(args.input_dir.glob('*.pt'))
    if not files:
        raise FileNotFoundError(f'No XFuse gene tensors in {args.input_dir}')
    for path in files:
        values = torch.load(path, map_location='cpu', weights_only=False)
        if isinstance(values, torch.Tensor):
            values = values.detach().cpu().numpy()
        values = np.asarray(values)
        if values.ndim == 3:
            values = values.mean(axis=0)
        elif values.ndim != 2:
            raise ValueError(f'Unexpected XFuse map dimensions: {values.shape}')
        with gzip.open(path.with_suffix('.pkl.gz'), 'wb', compresslevel=6) as handle:
            pickle.dump(values, handle)
    print(f'Exported {len(files)} XFuse gene maps; source tensors retained.')


if __name__ == '__main__':
    main()
