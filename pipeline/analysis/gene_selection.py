"""Select in-sample and cross-resolution genes and calculate overlap statistics."""
from __future__ import annotations
import argparse
import os
import glob
import re
from pathlib import Path
from typing import Dict, List, Tuple, Set
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import fisher_exact

DATASET_MAPPING = [('Xe_hsa_breast_c', 'Dataset 1'), ('Xe_hsa_colon_h', 'Dataset 2'), ('Xe_hsa_kidney_h', 'Dataset 3'), ('Xe_hsa_heart_h', 'Dataset 4'), ('Xe_hsa_skin_h', 'Dataset 5'), ('Xe_hsa_liver_h', 'Dataset 6'), ('ViHD_mmu_embryo_h', 'Dataset 7'), ('ViHD_hsa_pancreas_h', 'Dataset 8'), ('ViHD_mmu_brain_h', 'Dataset 9'), ('ViHD_hsa_colon_c', 'Dataset 10'), ('ViHD_mmu_kidney_h', 'Dataset 11'), ('cos_hsa_sarcoma_c', 'Dataset 12'), ('cos_hsa_prostate_c', 'Dataset 13')]
DATASET_NAME_TO_ID = dict(DATASET_MAPPING)
DATASET_ID_TO_NAME = {v: k for k, v in DATASET_MAPPING}
DATASETS_WITH_RESOLUTION = [f'Dataset {i}' for i in range(1, 12)]

def compute_parnorm(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate normalized mean metric ranks; lower scores indicate better performance."""
    df = df.copy()
    for col in ['PCC', 'SCC', 'KLD', 'JSD']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['SCC_rank'] = df['SCC'].rank(ascending=False, method='average', na_option='bottom')
    df['PCC_rank'] = df['PCC'].rank(ascending=False, method='average', na_option='bottom')
    df['KLD_rank'] = df['KLD'].rank(ascending=True, method='average', na_option='bottom')
    df['JSD_rank'] = df['JSD'].rank(ascending=True, method='average', na_option='bottom')
    rank_cols = ['SCC_rank', 'PCC_rank', 'KLD_rank', 'JSD_rank']
    df['PAR'] = df[rank_cols].mean(axis=1)
    scaler = MinMaxScaler()
    par_values = df[['PAR']].values
    if np.std(par_values) > 1e-10:
        df['PARnorm'] = scaler.fit_transform(par_values).flatten()
    else:
        df['PARnorm'] = 0.5
    return df

def determine_top_n(dataset_name: str, n_genes: int) -> int:
    """Determine the target gene count from dataset type and available gene coverage."""
    if dataset_name.startswith('Xe'):
        return max(30, min(60, int(0.12 * n_genes)))
    elif dataset_name.startswith('ViHD'):
        return max(80, min(120, int(0.1 * n_genes)))
    else:
        return max(30, min(80, int(0.1 * n_genes)))

def load_f2_data(base_dir: str) -> pd.DataFrame:
    """Load replicate-averaged in-sample metric tables."""
    pattern = re.compile('^(.+?)_(iSTAR|scstGCN|TESLA|XFuse|SuperST_t|SuperST_t_g)\\.csv$')
    all_data = []
    for filepath in glob.glob(os.path.join(base_dir, '*.csv')):
        filename = os.path.basename(filepath)
        match = pattern.match(filename)
        if not match:
            continue
        dataset_name = match.group(1)
        method = match.group(2)
        if dataset_name not in DATASET_NAME_TO_ID:
            continue
        df = pd.read_csv(filepath)
        if 'gene' not in df.columns:
            continue
        df = compute_parnorm(df)
        df['dataset_name'] = dataset_name
        df['dataset_id'] = DATASET_NAME_TO_ID[dataset_name]
        df['method'] = method
        all_data.append(df)
    return pd.concat(all_data, ignore_index=True)

def identify_top3_methods(f2_data: pd.DataFrame) -> Dict[str, List[str]]:
    """Select three method outputs using composite scores from PCC, SCC, JSD and KLD."""
    top3_methods = {}
    for dataset_id in f2_data['dataset_id'].unique():
        ds_data = f2_data[f2_data['dataset_id'] == dataset_id].copy()
        method_scores = []
        for method in ds_data['method'].unique():
            method_data = ds_data[ds_data['method'] == method]
            pcc_med = method_data['PCC'].median()
            scc_med = method_data['SCC'].median()
            jsd_med = method_data['JSD'].median()
            kld_med = method_data['KLD'].median()
            pcc_norm = (pcc_med - ds_data['PCC'].min()) / (ds_data['PCC'].max() - ds_data['PCC'].min() + 1e-10)
            scc_norm = (scc_med - ds_data['SCC'].min()) / (ds_data['SCC'].max() - ds_data['SCC'].min() + 1e-10)
            jsd_norm = 1 - (jsd_med - ds_data['JSD'].min()) / (ds_data['JSD'].max() - ds_data['JSD'].min() + 1e-10)
            kld_norm = 1 - (kld_med - ds_data['KLD'].min()) / (ds_data['KLD'].max() - ds_data['KLD'].min() + 1e-10)
            composite_score = (pcc_norm + scc_norm + jsd_norm + kld_norm) / 4.0
            method_scores.append({'method': method, 'score': composite_score, 'pcc': pcc_med, 'scc': scc_med, 'jsd': jsd_med, 'kld': kld_med})
        method_scores_df = pd.DataFrame(method_scores).sort_values('score', ascending=False)
        top3 = method_scores_df.head(3)['method'].tolist()
        top3_methods[dataset_id] = top3
    return top3_methods

def identify_top_performing_genes(f2_data: pd.DataFrame, top3_methods: Dict[str, List[str]], percentile_corr: float=0.4, percentile_div: float=0.6, min_pcc: float=0.05, min_scc: float=0.05, relax_factors: List[float]=None) -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    """
    Identify top-performing genes using only the top-3 methods per dataset.

    Procedure (per dataset):
      1) Aggregate gene-level median metrics across the top-3 methods:
         PCC_median, SCC_median, JSD_median, KLD_median, PARnorm_median, and n_methods.
      2) Candidate filtering (value-based):
         - PCC_median >= max(min_pcc, q_corr(PCC_median))
         - SCC_median >= max(min_scc, q_corr(SCC_median))
         - JSD_median <= q_div(JSD_median)
         - KLD_median <= q_div(KLD_median)
         where q_corr = percentile_corr (default 0.40), q_div = percentile_div (default 0.60).
      3) If candidates < top_n, gradually relax thresholds:
         - For PCC/SCC (higher better): threshold *= relax_factor
         - For JSD/KLD (lower better): threshold /= relax_factor
         - relax_factor = 0 means fully relaxed (no constraints).
      4) Rank candidates by PARnorm_median (ascending) and select top_n genes.

    Returns:
      - A dataframe with per-gene stats and selection flags per dataset
      - A dict mapping dataset_id -> list of selected genes
    """
    if relax_factors is None:
        relax_factors = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
    filtered_chunks = []
    for dataset_id, methods in top3_methods.items():
        ds = f2_data[(f2_data['dataset_id'] == dataset_id) & f2_data['method'].isin(methods)].copy()
        if not ds.empty:
            filtered_chunks.append(ds)
    if not filtered_chunks:
        return (pd.DataFrame(), {})
    f2_filtered = pd.concat(filtered_chunks, ignore_index=True)
    gene_stats = f2_filtered.groupby(['dataset_id', 'dataset_name', 'gene'], as_index=False).agg(n_methods=('method', 'nunique'), PCC_median=('PCC', 'median'), SCC_median=('SCC', 'median'), KLD_median=('KLD', 'median'), JSD_median=('JSD', 'median'), PARnorm_median=('PARnorm', 'median'))
    results = []
    selected_genes: Dict[str, List[str]] = {}
    for dataset_id in gene_stats['dataset_id'].unique():
        ds_data = gene_stats[gene_stats['dataset_id'] == dataset_id].copy()
        if ds_data.empty:
            selected_genes[dataset_id] = []
            continue
        dataset_name = ds_data['dataset_name'].iloc[0]
        n_genes = len(ds_data)
        top_n = determine_top_n(dataset_name, n_genes)
        pcc_thr = max(min_pcc, float(ds_data['PCC_median'].quantile(percentile_corr)))
        scc_thr = max(min_scc, float(ds_data['SCC_median'].quantile(percentile_corr)))
        jsd_thr = float(ds_data['JSD_median'].quantile(percentile_div))
        kld_thr = float(ds_data['KLD_median'].quantile(percentile_div))
        candidates = ds_data.copy()
        used_relax = None
        for rf in relax_factors:
            if rf > 0:
                adj_pcc = pcc_thr * rf
                adj_scc = scc_thr * rf
                adj_jsd = jsd_thr / rf
                adj_kld = kld_thr / rf
            else:
                adj_pcc = -np.inf
                adj_scc = -np.inf
                adj_jsd = np.inf
                adj_kld = np.inf
            cand = ds_data[(ds_data['PCC_median'] >= adj_pcc) & (ds_data['SCC_median'] >= adj_scc) & (ds_data['JSD_median'] <= adj_jsd) & (ds_data['KLD_median'] <= adj_kld)].copy()
            candidates = cand
            used_relax = rf
            if len(candidates) >= top_n:
                break
        candidates = candidates.sort_values('PARnorm_median', ascending=True)
        selected = candidates.head(top_n)['gene'].tolist()
        ds_data['selected_top_performing'] = ds_data['gene'].isin(selected)
        ds_data['top_n'] = top_n
        ds_data['rank'] = ds_data['PARnorm_median'].rank(method='min', ascending=True)
        ds_data['pcc_thr_base'] = pcc_thr
        ds_data['scc_thr_base'] = scc_thr
        ds_data['jsd_thr_base'] = jsd_thr
        ds_data['kld_thr_base'] = kld_thr
        ds_data['relax_factor_used'] = used_relax
        ds_data['n_candidates_final'] = len(candidates)
        results.append(ds_data)
        selected_genes[dataset_id] = selected
    return (pd.concat(results, ignore_index=True), selected_genes)

def load_f5_data(base_dir: str, method: str='iSTAR') -> pd.DataFrame:
    """Load resolution-dependent metric summaries."""
    all_data = []
    for dataset_name, dataset_id in DATASET_MAPPING:
        if dataset_id not in DATASETS_WITH_RESOLUTION:
            continue
        filepath = os.path.join(base_dir, method, dataset_name, 'res_mean.csv')
        if not os.path.exists(filepath):
            print(f'  Warning: File not found: {filepath}')
            continue
        try:
            df = pd.read_csv(filepath)
            required_cols = ['gene', 'simulation', 'resolution']
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                print(f'  Warning: Missing columns {missing_cols} in {filepath}')
                continue
            df['dataset_name'] = dataset_name
            df['dataset_id'] = dataset_id
            all_data.append(df)
            print(f'  Loaded {len(df)} records from {dataset_id} ({dataset_name})')
        except Exception as e:
            print(f'  Error loading {filepath}: {str(e)}')
            continue
    if len(all_data) == 0:
        raise ValueError(f'No data files found in {base_dir}/{method}/\nExpected structure: {base_dir}/{method}/{{dataset_name}}/res_mean.csv\nPlease check:\n  1. Base directory exists: {os.path.exists(base_dir)}\n  2. Method subdirectory exists: {os.path.exists(os.path.join(base_dir, method))}\n  3. Dataset subdirectories contain res_mean.csv files')
    return pd.concat(all_data, ignore_index=True)

def identify_cross_resolution_robust_genes(f5_data: pd.DataFrame, top_performing_counts: Dict[str, int], quality_percentile: float=20.0, stability_cv_max: float=0.5, min_conditions: int=6, count_ratio: float=1.0) -> Tuple[pd.DataFrame, Dict[str, List[str]], pd.DataFrame]:
    """Select genes with favorable performance and stability across resolutions.

    Filter by the quality percentile, CV threshold and minimum condition count.
    count_ratio sets the target count relative to the in-sample selected set.
    Return gene statistics, selected genes and condition-level observations.
    """
    print('\n=== Identifying Cross-Resolution Robust Genes ===')
    print(f'Quality criterion: Top {quality_percentile}% within each dataset')
    print(f'Stability criterion: CV < {stability_cv_max}')
    print(f'Count ratio: {count_ratio * 100}% of top-performing genes')
    print('\nStep 1: Computing PARnorm for all conditions...')
    condition_groups = []
    for (ds_id, sim, res), group in f5_data.groupby(['dataset_id', 'simulation', 'resolution']):
        group = compute_parnorm(group)
        condition_groups.append(group)
    f5_with_parnorm = pd.concat(condition_groups, ignore_index=True)
    print(f"  Total conditions: {len(f5_with_parnorm.groupby(['dataset_id', 'simulation', 'resolution']))}")
    print('\nStep 2: Aggregating gene-level statistics...')
    gene_stats = f5_with_parnorm.groupby(['dataset_id', 'dataset_name', 'gene']).agg({'PARnorm': ['mean', 'std', lambda x: np.percentile(x.dropna(), 75), 'count'], 'PCC': 'median', 'SCC': 'median'}).reset_index()
    gene_stats.columns = ['dataset_id', 'dataset_name', 'gene', 'PARnorm_mean', 'PARnorm_std', 'PARnorm_q75', 'n_conditions', 'PCC_median', 'SCC_median']
    gene_stats['CV'] = gene_stats['PARnorm_std'] / (gene_stats['PARnorm_mean'] + 1e-06)
    print('\nStep 3: Selecting genes for each dataset...')
    results = []
    selected_genes = {}
    for dataset_id in DATASETS_WITH_RESOLUTION:
        ds_data = gene_stats[gene_stats['dataset_id'] == dataset_id].copy()
        if ds_data.empty:
            print(f'  {dataset_id}: No data available')
            selected_genes[dataset_id] = []
            continue
        dataset_name = ds_data['dataset_name'].iloc[0]
        quality_threshold = np.percentile(ds_data['PARnorm_mean'], quality_percentile)
        quality_pass = ds_data['PARnorm_mean'] <= quality_threshold
        stability_pass = ds_data['CV'] <= stability_cv_max
        sufficient_data = ds_data['n_conditions'] >= min_conditions
        ds_data['QualityScore'] = 1 - ds_data['PARnorm_mean']
        ds_data['StabilityScore'] = 1 - ds_data['CV'].clip(0, 1)
        ds_data['CompositeScore'] = 0.7 * ds_data['QualityScore'] + 0.3 * ds_data['StabilityScore']
        is_candidate = quality_pass & stability_pass & sufficient_data
        candidates = ds_data[is_candidate].copy()
        if dataset_id in top_performing_counts:
            target_n = int(top_performing_counts[dataset_id] * count_ratio)
        else:
            target_n = determine_top_n(dataset_name, len(ds_data))
        actual_n = min(target_n, len(candidates))
        selected = candidates.sort_values('CompositeScore', ascending=False).head(actual_n)['gene'].tolist()
        ds_data['selected_cross_resolution'] = ds_data['gene'].isin(selected)
        ds_data['quality_threshold'] = quality_threshold
        ds_data['target_n'] = target_n
        ds_data['actual_n'] = actual_n
        ds_data['n_candidates'] = len(candidates)
        print(f'  {dataset_id}:')
        print(f'    Total genes: {len(ds_data)}')
        print(f'    Quality threshold (PARnorm_mean ≤ {quality_threshold:.3f}): {quality_pass.sum()} genes')
        print(f'    Stability pass (CV ≤ {stability_cv_max}): {stability_pass.sum()} genes')
        print(f'    Sufficient data (≥ {min_conditions} conditions): {sufficient_data.sum()} genes')
        print(f'    Candidates (all criteria): {len(candidates)} genes')
        print(f'    Target count: {target_n}, Selected: {actual_n} genes')
        if actual_n < target_n:
            print(f'    ⚠ Warning: Only {actual_n}/{target_n} genes meet all criteria')
        results.append(ds_data)
        selected_genes[dataset_id] = selected
    final_stats = pd.concat(results, ignore_index=True)
    print('\n=== Summary ===')
    total_selected = sum((len(genes) for genes in selected_genes.values()))
    print(f'Total cross-resolution robust genes selected: {total_selected}')
    print(f'Average per dataset: {total_selected / len(DATASETS_WITH_RESOLUTION):.1f}')
    return (final_stats, selected_genes, f5_with_parnorm)

def compute_overlap_stats(top_genes: Dict[str, List[str]], cross_genes: Dict[str, List[str]], all_genes_count: Dict[str, int]) -> pd.DataFrame:
    """Calculate gene-set overlap statistics for each dataset."""
    stats = []
    for dataset_id in DATASETS_WITH_RESOLUTION:
        set_a = set(top_genes.get(dataset_id, []))
        set_b = set(cross_genes.get(dataset_id, []))
        overlap = set_a & set_b
        union = set_a | set_b
        jaccard = len(overlap) / len(union) if len(union) > 0 else 0
        n_total = all_genes_count.get(dataset_id, 1000)
        a = len(overlap)
        b = len(set_a - set_b)
        c = len(set_b - set_a)
        d = max(0, n_total - a - b - c)
        try:
            odds, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        except:
            odds, pval = (np.nan, np.nan)
        stats.append({'dataset_id': dataset_id, 'n_top_performing': len(set_a), 'n_cross_resolution': len(set_b), 'n_overlap': len(overlap), 'jaccard': jaccard, 'overlap_pct_top': len(overlap) / len(set_a) * 100 if len(set_a) > 0 else 0, 'overlap_pct_cross': len(overlap) / len(set_b) * 100 if len(set_b) > 0 else 0, 'fisher_odds': odds, 'fisher_p': pval})
    return pd.DataFrame(stats)

def export_supplementary_tables(top_stats: pd.DataFrame, cross_stats: pd.DataFrame, overlap_stats: pd.DataFrame, top3_methods: Dict[str, List[str]], output_dir: str, quality_percentile: float=20.0, stability_cv_max: float=0.5, count_ratio: float=1.0):
    """Export gene statistics, selection labels and overlap results as workbooks."""
    os.makedirs(output_dir, exist_ok=True)
    table_s1 = top_stats[['dataset_id', 'dataset_name', 'gene', 'n_methods', 'PCC_median', 'SCC_median', 'KLD_median', 'JSD_median', 'PARnorm_median', 'rank', 'selected_top_performing']].copy()
    table_s1 = table_s1.sort_values(['dataset_id', 'selected_top_performing', 'PARnorm_median'], ascending=[True, False, True])
    col_desc_s1 = pd.DataFrame({'Column': table_s1.columns, 'Description': ['Dataset identifier (Dataset1-Dataset13)', 'Original dataset name', 'Gene symbol', 'Number of top-3 methods that evaluated this gene', 'Median PCC across top-3 methods (higher is better)', 'Median SCC across top-3 methods (higher is better)', 'Median KLD across top-3 methods (lower is better)', 'Median JSD across top-3 methods (lower is better)', 'Median PARnorm across top-3 methods (0-1 scale, lower is better)', 'Gene rank within dataset based on PARnorm_median', 'Whether selected as top-performing gene (TRUE/FALSE)']})
    top3_info = pd.DataFrame([{'dataset_id': ds_id, 'top3_methods': ', '.join(methods)} for ds_id, methods in top3_methods.items()])
    with pd.ExcelWriter(os.path.join(output_dir, 'Gu_Supplementary_Table_3.xlsx'), engine='openpyxl') as writer:
        col_desc_s1.to_excel(writer, sheet_name='Column_Descriptions', index=False)
        top3_info.to_excel(writer, sheet_name='Top3_Methods_Per_Dataset', index=False)
        table_s1.to_excel(writer, sheet_name='Gene_Statistics', index=False)
    print(f"Table 3 exported: {len(table_s1[table_s1['selected_top_performing']])} top-performing genes")
    table_s2 = cross_stats[['dataset_id', 'dataset_name', 'gene', 'n_conditions', 'PARnorm_mean', 'PARnorm_std', 'PARnorm_q75', 'CV', 'PCC_median', 'SCC_median', 'QualityScore', 'StabilityScore', 'CompositeScore', 'quality_threshold', 'n_candidates', 'selected_cross_resolution']].copy()
    table_s2 = table_s2.sort_values(['dataset_id', 'selected_cross_resolution', 'CompositeScore'], ascending=[True, False, False])
    col_desc_s2 = pd.DataFrame({'Column': table_s2.columns, 'Description': ['Dataset identifier (Dataset1-Dataset11, only datasets with resolution variations)', 'Original dataset name', 'Gene symbol', 'Number of (simulation, resolution) conditions evaluated', 'Mean PARnorm across all conditions (lower indicates better average performance)', 'Standard deviation of PARnorm across conditions (lower indicates more stable performance)', '75th percentile of PARnorm (upper-quartile relative performance, lower is better)', f'Coefficient of variation (CV = std/mean, lower indicates more stable performance; threshold: {stability_cv_max})', 'Median Pearson correlation across conditions', 'Median Spearman correlation across conditions', 'Quality score: 1 - PARnorm_mean (higher is better, range 0-1)', 'Stability score: 1 - clip(CV, 0, 1) (higher is better, range 0-1)', 'Composite score: 0.7×QualityScore + 0.3×StabilityScore (higher is better, range 0-1)', f'Quality threshold for this dataset (PARnorm_mean at {quality_percentile}th percentile)', 'Number of candidate genes meeting all filtering criteria', 'Whether selected as cross-resolution robust gene (TRUE/FALSE)']})
    selection_criteria = pd.DataFrame({'Criterion': ['Quality criterion', 'Stability criterion', 'Data sufficiency', 'Count control', 'Ranking method'], 'Description': [f'PARnorm_mean ≤ {quality_percentile}th percentile within each dataset', f'Coefficient of variation (CV) ≤ {stability_cv_max}', 'At least 6 conditions evaluated', f'Select approximately {count_ratio * 100}% of top-performing gene count', 'Rank by composite score (70% quality + 30% stability)']})
    with pd.ExcelWriter(os.path.join(output_dir, 'Gu_Supplementary_Table_4.xlsx'), engine='openpyxl') as writer:
        col_desc_s2.to_excel(writer, sheet_name='Column_Descriptions', index=False)
        selection_criteria.to_excel(writer, sheet_name='Selection_Criteria', index=False)
        table_s2.to_excel(writer, sheet_name='Gene_Statistics', index=False)
    print(f"Table 4 exported: {len(table_s2[table_s2['selected_cross_resolution']])} cross-resolution robust genes")
    col_desc_s3 = pd.DataFrame({'Column': overlap_stats.columns, 'Description': ['Dataset identifier (Dataset1-Dataset11)', 'Number of genes selected as top-performing', 'Number of genes selected as cross-resolution robust', 'Number of genes in both sets (overlap)', 'Jaccard index (overlap size / union size, measures similarity between two sets, range 0-1)', 'Overlap as percentage of top-performing set', 'Overlap as percentage of cross-resolution set', "Fisher's exact test odds ratio (measures enrichment of overlap, >1 indicates positive association)", "Fisher's exact test p-value (one-sided test for enrichment, p<0.05 indicates significant overlap)"]})
    statistical_notes = pd.DataFrame({'Metric': ['Jaccard Index', 'Overlap Percentage', "Fisher's Exact Test", 'Significance Levels'], 'Interpretation': ['Measures similarity between two sets. J=0: no overlap; J=1: identical sets; J>0.5: substantial overlap', 'Shows what proportion of each gene set is shared. High percentages indicate strong consistency', 'Tests whether overlap is greater than expected by chance. Odds ratio >1 with p<0.05 indicates significant enrichment', '*** p<0.001 (highly significant); ** p<0.01 (very significant); * p<0.05 (significant)']})
    with pd.ExcelWriter(os.path.join(output_dir, 'Gu_Supplementary_Table_5.xlsx'), engine='openpyxl') as writer:
        col_desc_s3.to_excel(writer, sheet_name='Column_Descriptions', index=False)
        statistical_notes.to_excel(writer, sheet_name='Statistical_Notes', index=False)
        overlap_stats.to_excel(writer, sheet_name='Overlap_Statistics', index=False)
    print(f'Table 5 exported: overlap analysis for {len(overlap_stats)} datasets')
    method_summary = pd.DataFrame({'Analysis': ['Top-performing genes', 'Cross-resolution robust genes', 'Overlap analysis'], 'Method': ['Select genes with lowest median PARnorm across top-3 methods', f'Three-step filtering: (1) Quality: top {quality_percentile}% PARnorm_mean; (2) Stability: CV≤{stability_cv_max}; (3) Count: ~{count_ratio * 100}% of top-performing count', "Compare two gene sets using Jaccard index and Fisher's exact test"], 'Output_Table': ['Gu_Supplementary_Table_3.xlsx', 'Gu_Supplementary_Table_4.xlsx', 'Gu_Supplementary_Table_5.xlsx']})
    method_summary.to_excel(os.path.join(output_dir, 'README_SupplementaryTables.xlsx'), sheet_name='Method_Summary', index=False)
    print(f'Method summary exported: README_SupplementaryTables.xlsx')
    print(f'\nAll supplementary tables exported to: {output_dir}')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--insample-dir',type=Path,required=True,help='Replicate-averaged dataset_method.csv files')
    p.add_argument('--resolution-dir',type=Path,required=True,help='iSTAR/dataset_name/res_mean.csv directories')
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args()
    f2=load_f2_data(str(args.insample_dir))
    expected=set(DATASET_NAME_TO_ID.values())
    if set(f2.dataset_id)!=expected:raise ValueError('In-sample inputs must cover all 13 study datasets')
    for dataset, group in f2.groupby('dataset_id'):
        if set(group.method)!={'XFuse','TESLA','iSTAR','scstGCN','SuperST_t','SuperST_t_g'}:
            raise ValueError(f'Missing method outputs for {dataset}')
    if f2.duplicated(['dataset_id','method','gene']).any():raise ValueError('Duplicate in-sample dataset-method-gene rows')
    top3=identify_top3_methods(f2)
    top_stats,top_genes=identify_top_performing_genes(f2,top3)
    f5=load_f5_data(str(args.resolution_dir),method='iSTAR')
    if set(f5.dataset_id)!=set(DATASETS_WITH_RESOLUTION):raise ValueError('Resolution inputs must cover all 11 evaluable study datasets')
    if f5.duplicated(['dataset_id','simulation','resolution','gene']).any():raise ValueError('Duplicate resolution-condition-gene rows')
    cross_stats,cross_genes,conditions=identify_cross_resolution_robust_genes(f5,{d:len(g) for d,g in top_genes.items()})
    overlap=compute_overlap_stats(top_genes,cross_genes,{d:len(g) for d,g in top_stats.groupby('dataset_id')})
    args.output_dir.mkdir(parents=True,exist_ok=True)
    export_supplementary_tables(top_stats,cross_stats,overlap,top3,str(args.output_dir))
    top_stats.to_csv(args.output_dir/'top_gene_selection_audit.csv',index=False)
    cross_stats.to_csv(args.output_dir/'cross_resolution_selection_audit.csv',index=False)
    print('Saved final Tables 3, 4 and 5; numerical selection rules are unchanged from the supplied source.')

if __name__=='__main__':main()
