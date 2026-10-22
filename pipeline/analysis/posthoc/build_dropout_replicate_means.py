"""Aggregate per-gene dropout metrics into replicate-level means.

The aggregation step preserves the original per-gene logic;
rep1/rep2/rep3 correspond to seeds 123/124/125.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


OUT = Path(__file__).resolve().parent
DATASETS = ['Xe_hsa_breast_c','Xe_hsa_colon_h','Xe_hsa_kidney_h',
            'Xe_hsa_heart_h','Xe_hsa_skin_h','Xe_hsa_liver_h']
METHODS = ['XFuse','TESLA','iSTAR','scstGCN','SuperST_t','SuperST_t_g']
RATIOS = [0.1,0.2,0.3,0.4,0.5]
METRICS = ['PCC','SCC','JSD','KLD']


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input-root',type=Path,required=True,help='Metric CSV root: method/dataset/drop_out/drop_out_0.1/rep1.csv')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    rows=[]
    for dataset in DATASETS:
        for method in METHODS:
            if method=='scstGCN' and dataset=='Xe_hsa_skin_h':
                continue  # All 15 attempted outputs were degenerate and non-evaluable; no metric values are imputed.
            for ratio in RATIOS:
                for rep in (1,2,3):
                    path=args.input_root/method/dataset/'drop_out'/f'drop_out_{ratio:.1f}'/f'rep{rep}.csv'
                    if not path.exists():
                        raise FileNotFoundError(path)
                    df=pd.read_csv(path,usecols=METRICS)
                    if df.empty:
                        raise ValueError(f'Empty metric table: {path}')
                    for metric in METRICS:
                        vals=pd.to_numeric(df[metric],errors='coerce').replace([np.inf,-np.inf],np.nan)
                        if not vals.notna().any():
                            raise ValueError(f'No finite {metric} values: {path}')
                        rows.append(dict(dataset=dataset,method=method,metric=metric,
                            dropout_ratio=ratio,replicate=rep,
                            n_valid_genes=int(vals.notna().sum()),
                            replicate_gene_mean=float(vals.mean())))
    replicate=pd.DataFrame(rows)
    assert len(replicate)==2100
    args.output.parent.mkdir(parents=True,exist_ok=True)
    replicate.to_csv(args.output,index=False,float_format='%.17g')
    print('Saved',args.output,'with',len(replicate),'rows')

if __name__=='__main__':
    main()
