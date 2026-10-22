"""Command-line entry points for the reported SRBench analyses.

Run from the archive root: python -m pipeline.analysis.cli --help
The numerical implementations live in the companion modules.
"""
from __future__ import annotations
import argparse
import importlib
import json
from pathlib import Path

METRICS = ['PCC', 'SCC', 'JSD', 'KLD', 'SSIM', 'PSNR']
METHODS = ['XFuse', 'TESLA', 'iSTAR', 'scstGCN', 'SuperST_t', 'SuperST_t_g']


def load(name):
    return importlib.import_module('pipeline.analysis.' + name)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='task', required=True)
    q = sub.add_parser('export-source', help='Export existing replicate-averaged workbook rows to analysis CSVs')
    q.add_argument('--workbook', type=Path, required=True)
    q.add_argument('--output-dir', type=Path, required=True)
    q.add_argument('--pipeline-gene-ids', action='store_true', help='Apply the documented CosMx slash-to-underscore filename mapping and record it')
    q = sub.add_parser('insample', help='Aggregate three runs and compute registered robustness')
    q.add_argument('--metric-root', type=Path, required=True)
    q.add_argument('--config', type=Path, required=True, help='JSON: method_datasets mapping')
    q.add_argument('--output-dir', type=Path, required=True)
    q = sub.add_parser('rank-insample', help='Rank methods within each dataset')
    q.add_argument('--input-dir', type=Path, required=True)
    q.add_argument('--datasets', nargs='+', required=True)
    q.add_argument('--output-dir', type=Path, required=True)
    q.add_argument('--plot', action='store_true')
    q = sub.add_parser('hvg', help='HVG/non-HVG figures and two-sided Mann–Whitney U-test results')
    q.add_argument('--input-dir', type=Path, required=True)
    q.add_argument('--dataset', required=True)
    q.add_argument('--methods', nargs='+', default=METHODS)
    q.add_argument('--output-dir', type=Path, required=True)
    q = sub.add_parser('strata', help='Expression-by-Moran-I heatmaps, means and precise gene counts')
    q.add_argument('--input-dir', type=Path, required=True)
    q.add_argument('--dataset', required=True)
    q.add_argument('--methods', nargs='+', default=METHODS)
    q.add_argument('--output-dir', type=Path, required=True)
    q = sub.add_parser('gene-effects', help='Fig. 3d standardized high-minus-low contrasts')
    q.add_argument('--input-dir', type=Path, required=True)
    q.add_argument('--output-dir', type=Path, required=True)
    q = sub.add_parser('top-genes', help='In-sample top-gene selection without requiring resolution inputs')
    q.add_argument('--input-dir', type=Path, required=True)
    q.add_argument('--output-dir', type=Path, required=True)
    q = sub.add_parser('fov', help='FOV rankings, optional location maps and size box plots')
    q.add_argument('--config', type=Path, required=True, help='fov_dirs, method, optional he_image')
    q.add_argument('--output-dir', type=Path, required=True)
    q = sub.add_parser('adjacent', help='Gene-paired within-slide/cross-slide figures and Wilcoxon tests')
    q.add_argument('--input-dir', type=Path, required=True)
    q.add_argument('--methods', nargs='+', required=True)
    q.add_argument('--output-dir', type=Path, required=True)
    q = sub.add_parser('resolution', help='Replicate means, ranks and mean-only resolution curves')
    q.add_argument('--input-dir', type=Path, required=True, help='rep*/scheme/parameter_value.csv')
    q.add_argument('--method', required=True, choices=METHODS)
    q.add_argument('--dataset', required=True)
    q.add_argument('--output-dir', type=Path, required=True)
    q = sub.add_parser('usability', help='Normalize actual usability observations and draw Fig. 6')
    q.add_argument('--input', type=Path, required=True)
    q.add_argument('--output-dir', type=Path, required=True)
    args = p.parse_args(argv)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    if args.task == 'export-source':
        load('source_tables').export_workbook(args.workbook, out, pipeline_gene_ids=args.pipeline_gene_ids)
        return
    import numpy as np
    import pandas as pd
    if args.task == 'insample':
        config = json.loads(args.config.read_text(encoding='utf-8'))
        load('fig2_insample').merge_by_method_datasets_and_compute_robustness(
            str(args.metric_root), str(out), config['method_datasets'])
    elif args.task == 'rank-insample':
        module = load('fig2_insample')
        df = module.build_df_for_plot(str(args.input_dir), args.datasets, methods=METHODS)
        df.to_csv(out / 'insample_gene_metrics_and_ranks.csv', index=False, float_format='%.17g')
        df[['dataset', 'method', 'PArank', 'SArank']].drop_duplicates().to_csv(out / 'method_ranks.csv', index=False)
        if args.plot:
            module.plot_panel_from_df(df, method_order=METHODS, plot_name=str(out / 'insample_overview.pdf'), show=False)
    elif args.task in {'hvg', 'strata'}:
        from scipy.stats import mannwhitneyu
        rows = []
        levels = ['Low', 'Medium', 'High']
        for method in args.methods:
            path = args.input_dir / f'{args.dataset}_{method}_classified.csv'
            df = pd.read_csv(path)
            if df['gene'].duplicated().any():
                raise ValueError(f'Duplicate gene identifiers: {path}')
            if args.task == 'hvg':
                for metric in METRICS:
                    a = pd.to_numeric(df.loc[df.hvg_status == 'HVG', metric], errors='coerce').dropna()
                    b = pd.to_numeric(df.loc[df.hvg_status == 'non-HVG', metric], errors='coerce').dropna()
                    if not np.isfinite(a).all() or not np.isfinite(b).all():
                        raise ValueError(f'Non-finite {metric} values: {path}')
                    statistic, pv = mannwhitneyu(a, b, alternative='two-sided') if len(a) and len(b) else (np.nan, np.nan)
                    rows.append(dict(dataset=args.dataset, method=method, metric=metric,
                                     n_HVG=len(a), n_non_HVG=len(b), U=statistic, P_value=pv,
                                     test='two-sided Mann-Whitney U', adjustment='none'))
                load('fig3_hvg').plot_hvg_comparison(method, args.dataset, str(args.input_dir), str(out),
                    box_width=.4, subplot_spacing=.07, spine_linewidth=.4, box_linewidth=.4)
            else:
                for metric in METRICS:
                    load('fig3_heatmaps').plot_3x3_metric_heatmap(method, args.dataset, metric, str(args.input_dir), str(out))
                    for spatial in levels:
                        for expression in levels:
                            values = pd.to_numeric(df.loc[
                                (df.Spatial_Variability_Level == spatial) &
                                (df.Gene_Expression_Level == expression), metric], errors='coerce')
                            values = values[np.isfinite(values)]
                            rows.append(dict(dataset=args.dataset, method=method, metric=metric,
                                             Spatial_Variability_Level=spatial,
                                             Gene_Expression_Level=expression,
                                             n_genes=len(values), mean_value=values.mean()))
        pd.DataFrame(rows).to_csv(out / f'{args.task}_statistics.csv', index=False, float_format='%.17g')
    elif args.task == 'gene-effects':
        load('fig3_summary').run_overall_summary(str(args.input_dir),
            str(out / 'gene_property_effects.csv'), str(out / 'gene_property_effects.pdf'))
    elif args.task == 'top-genes':
        module = load('gene_selection')
        frame = module.load_f2_data(str(args.input_dir))
        if set(frame.dataset_id) != set(module.DATASET_NAME_TO_ID.values()):
            raise ValueError('The study top-gene analysis requires all 13 datasets')
        if frame.duplicated(['dataset_id', 'method', 'gene']).any():
            raise ValueError('Duplicate dataset-method-gene rows')
        for dataset, group in frame.groupby('dataset_id'):
            if set(group.method) != set(METHODS):
                raise ValueError(f'Missing method outputs for {dataset}')
        top3 = module.identify_top3_methods(frame)
        statistics, selected = module.identify_top_performing_genes(frame, top3)
        statistics.to_csv(out / 'top_gene_statistics.csv', index=False, float_format='%.17g')
        (out / 'selected_top_genes.json').write_text(json.dumps(selected, indent=2) + '\n')
        (out / 'top3_method_outputs.json').write_text(json.dumps(top3, indent=2) + '\n')
    elif args.task == 'fov':
        config = json.loads(args.config.read_text(encoding='utf-8'))
        base = args.config.resolve().parent
        resolve = lambda x: str((base / x).resolve()) if not Path(x).is_absolute() else x
        dirs = {str(k): resolve(v) for k, v in config['fov_dirs'].items()}
        big_csv = out / 'fov_gene_metrics.csv'
        load('fig4_fov').prepare_rank_par_norm_fovsize_level(dirs, str(big_csv))
        boxes = load('fig4_fov_boxplots')
        frame, score = boxes.preprocess_PARnorm(str(big_csv))
        frame.to_csv(big_csv, index=False, float_format='%.17g')
        sizes = tuple(int(str(k).lower().replace('fov', '').replace('_', '')) for k in dirs)
        boxes.plot_box_violin_fov(config['method'], str(big_csv), [score], str(out), fov_size_list=sizes)
        if config.get('he_image'):
            from PIL import Image
            image = np.asarray(Image.open(resolve(config['he_image'])))
            load('fig4_fov').plot_fovs(image, str(big_csv), 'PARnorm', fig_name=str(out / 'fov_location_maps.pdf'))
    elif args.task == 'adjacent':
        from scipy.stats import wilcoxon
        rows = []
        for method in args.methods:
            a = pd.read_csv(args.input_dir / method / 'within_slice.csv')
            b = pd.read_csv(args.input_dir / method / 'cross_slice.csv')
            paired = a.merge(b, on='gene', suffixes=('_within', '_cross'), validate='one_to_one')
            for metric in METRICS:
                pair = paired[[metric + '_within', metric + '_cross']].replace([np.inf, -np.inf], np.nan).dropna()
                statistic, pv = wilcoxon(pair.iloc[:, 0], pair.iloc[:, 1], alternative='two-sided') if len(pair) else (np.nan, np.nan)
                rows.append(dict(method=method, metric=metric, n_pairs=len(pair), W=statistic,
                                 P_value=pv, test='two-sided Wilcoxon signed-rank', adjustment='none'))
            load('fig4_adjacent').plot_slice_comparison(method, str(args.input_dir), str(out),
                box_width=.4, subplot_spacing=.07, spine_linewidth=.4, box_linewidth=.4)
        pd.DataFrame(rows).to_csv(out / 'adjacent_statistics.csv', index=False, float_format='%.17g')
    elif args.task == 'resolution':
        files = sorted(args.input_dir.glob('rep*/*/*.csv'))
        if not files:
            raise FileNotFoundError(f'No rep*/scheme/*.csv files in {args.input_dir}')
        module = load('fig5_resolution')
        ranked = module.process_data(files, str(out), args.method, args.dataset)
        colors = dict(proportional='#0072B2', dense='#D55E00', fixed_spotDistance='#009E73', fixed_spotSize='#CC79A7')
        module.plot_all_metrics(args.dataset, res_rank=ranked, metrics=['PARnorm'] + METRICS,
                                method=args.method, plt_path=str(out), color_map=colors,
                                resolution_types=('spotdistance',) if args.method == 'TESLA' else ('spotsize', 'spotdistance'))
    elif args.task == 'usability':
        module = load('fig6_usability')
        raw = module._read_usability_table(str(args.input))
        module.compute_normalized_scores(raw).to_csv(out / 'normalized_usability.csv', float_format='%.17g')
        module.plot_usability_heatmap(raw, save_pdf=str(out / 'usability.pdf'), show=False)


if __name__ == '__main__':
    main()
