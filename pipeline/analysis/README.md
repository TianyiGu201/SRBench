# Running the analyses

Run commands from the archive root. Paths named `DATA` refer to supplied or regenerated study data; derived files are written under `RESULTS`. All examples use the Supplementary Table numbering in the manuscript.

## 1. Choose the input route

**Raw metric CSVs:** evaluate each prediction with `pipeline/metrics`, then use `insample` to average the three computational runs and calculate robustness.

```bash
python -m pipeline.analysis.cli insample --metric-root DATA/metric_csv --config configs/insample.example.json --output-dir RESULTS/merged
python -m pipeline.analysis.cli rank-insample --input-dir RESULTS/merged --datasets Xe_hsa_breast_c --output-dir RESULTS/ranks --plot
```

**Published gene-level Source Data:** export the existing `*_raw` workbook sheets. These values have already been averaged across runs. No extra averaging is done.

```bash
python -m pipeline.analysis.cli export-source --workbook DATA/Gu_Supplementary_Data_SuppFig6_95.xlsx --output-dir RESULTS/source_tables
```

This creates `merged/<dataset>_<method>.csv` and, where annotations are present, `classified/<dataset>_<method>_classified.csv`, plus a sheet-to-file manifest. Existing HVG labels and expression/Moran-I strata are retained. The exporter does not silently reclassify boundary genes or fill absent genes.

For the pipeline's file-safe CosMx gene identifiers, add `--pipeline-gene-ids`. This applies the existing preprocessing replacement `/` → `_`, retains the original target label in `source_gene`, and writes `gene_identifier_mapping.csv`. Use this option when generating Table 3 with the pipeline identifiers:

```bash
python -m pipeline.analysis.cli export-source --workbook DATA/Gu_Supplementary_Data_SuppFig6_95.xlsx --output-dir RESULTS/source_tables --pipeline-gene-ids
```

It does not change metric values or gene-group annotations.

The six output names are `XFuse`, `TESLA`, `iSTAR`, `scstGCN`, `SuperST_t` and `SuperST_t_g`. Exact dataset identifiers are recorded in `gene_selection.DATASET_NAME_TO_ID`.

## 2. Gene-level effects

Classified CSVs need `gene`, `dataset`, `method`, the six metrics, `hvg_status`, `morans_I`, `gt_counts`, `Gene_Expression_Level` and `Spatial_Variability_Level`.

```bash
python -m pipeline.analysis.cli hvg --input-dir RESULTS/source_tables/classified --dataset Xe_hsa_breast_c --methods iSTAR scstGCN --output-dir RESULTS/hvg
python -m pipeline.analysis.cli strata --input-dir RESULTS/source_tables/classified --dataset cos_hsa_sarcoma_c --methods XFuse TESLA iSTAR scstGCN SuperST_t SuperST_t_g --output-dir RESULTS/strata
python -m pipeline.analysis.cli gene-effects --input-dir RESULTS/source_tables/classified --output-dir RESULTS/gene_effects
```

- `hvg`: six-metric PDFs and two-sided Mann–Whitney U statistics, P values and per-group n in `hvg_statistics.csv`; these P values are unadjusted.
- `strata`: six metric heatmaps per method and `strata_statistics.csv` with all nine strata, including n=0 cells, finite gene counts and means.
- `gene-effects`: within-dataset/method z-score composites, high-minus-low marginal contrasts, `gene_property_effects.csv` and a Fig. 3d PDF. Missing strata are omitted separately in each marginal mean; missing contrasts remain missing. Large markers show across-dataset means and coordinate-wise IQRs.

`gene_characterization.py` computes annotations from platform inputs when the Scanpy/Squidpy environment is available. `gene_groups.classify_genes` applies the registered independent strata. The published workbook route preserves the existing labels instead of recomputing them.

## 3. FOV and adjacent-section analyses

FOV inputs are per-position gene-metric CSVs. Their filenames must contain the coordinate token `minx..._maxx..._miny..._maxy...`. Paths in the FOV configuration are relative to that JSON file, not to the current shell directory.

```bash
python -m pipeline.analysis.cli fov --config configs/fov.example.json --output-dir RESULTS/fov
python -m pipeline.analysis.cli adjacent --input-dir DATA/adjacent --methods iSTAR scstGCN --output-dir RESULTS/adjacent
```

FOV ranking has two different scopes, retained from the supplied study code: location maps rank pooled gene–FOV observations **within each FOV size**; the size-comparison box plot uses `PARnorm_all_fov`, ranked **across sizes**. The derived CSV is written to the output directory, not back into source inputs. Whole-FOV failures retain their rows, but their final PAR/PARnorm fields are blank and they do not contribute to plots or summaries. The existing ranking operations are retained before applying this status mask, preserving evaluable scores. An undefined individual metric is distinguished from a whole-FOV failure.

Adjacent inputs are `<method>/within_slice.csv` and `<method>/cross_slice.csv`. Genes are matched by identifier, never by row order. Outputs include a paired comparison PDF and `adjacent_statistics.csv` containing n, W and two-sided Wilcoxon signed-rank P values. Both test workflows use SciPy's automatic method selection. These P values are unadjusted. Gene counts do not represent independent tissue replicates.

## 4. Initial resolution

```bash
python -m pipeline.analysis.cli resolution --input-dir DATA/resolution/iSTAR/Xe_hsa_breast_c --method iSTAR --dataset Xe_hsa_breast_c --output-dir RESULTS/resolution/iSTAR/Xe_hsa_breast_c
```

Input layout is `rep1|rep2|rep3/<scheme>/spotSize_<value>.csv` or `.../fixed_spotSize/spotDistance_<value>.csv`. Each table includes `gene`, `gt_counts` and the six raw metrics. Supply all four evaluated schemes for non-TESLA outputs; TESLA uses only `fixed_spotSize`.

Outputs are `res_all.csv`, `res_mean.csv`, `res_rank.csv` and mean-only curves. The code averages computational runs before ranking gene/configuration rows within each scheme. The curves do **not** estimate or draw error bars or CIs.

## 5. Gene selection and Supplementary Tables 3–5

In-sample selection can be run independently of the resolution analyses:

```bash
python -m pipeline.analysis.cli top-genes --input-dir RESULTS/source_tables/merged --output-dir RESULTS/top_genes
```

This writes gene statistics, selected gene identifiers and selected method outputs. To generate all three formatted workbooks, supply both input collections:

```bash
python -m pipeline.analysis.gene_selection --insample-dir RESULTS/source_tables/merged --resolution-dir RESULTS/resolution --output-dir RESULTS/gene_selection
```

The full-study command requires all 13 datasets and six outputs for in-sample selection, plus `iSTAR/<dataset>/res_mean.csv` for the 11 resolution datasets. Provide the complete input collection for this command.

The supplied calculation rules are retained: composite scores select three method outputs; their per-gene medians undergo joint metric filtering with the actual dataset-specific threshold-relaxation sequence; candidates are ordered by median PARnorm. Cross-resolution selection ranks genes **within each configuration**, then summarizes ranks over configurations, applies quality, CV and condition-count criteria, and uses the retained composite ranking. This ranking scope differs from the resolution-curve ranking scope above. Overlap uses a one-sided Fisher exact test. The export contains Tables 3–5 and gene-level selection audits.

## 6. Post hoc absolute in-sample summaries: Table 2

```bash
python -m pipeline.analysis.posthoc.build_insample_absolute_summaries --fig2-source DATA/Gu_SourceData_Fig2.xlsx --suppfig5-source DATA/Gu_Supplementary_Data_SuppFig5.xlsx --output RESULTS/Gu_Supplementary_Table_2.xlsx
```

The builder reads the original source-workbook tables directly. It retains zero values, averages no extra replicate layer, reports finite-gene denominators, and measures coverage against each dataset's gene union. Mean squared PCC is the mean of squared gene-level PCC values, not the residual-based R-squared.

## 7. Registered dropout simulations and post hoc Table 6

```bash
python -m pipeline.analysis.posthoc.build_dropout_replicate_means --input-root DATA/metric_csv --output RESULTS/dropout/dropout_replicate_means_all_Xenium.csv
python -m pipeline.analysis.posthoc.build_dropout_paired_statistics --input RESULTS/dropout/dropout_replicate_means_all_Xenium.csv --output-dir RESULTS/dropout
python -m pipeline.analysis.posthoc.build_dropout_supplement_a4 --input-dir RESULTS/dropout --output-dir RESULTS/dropout_figures --mode raw
```

Metric inputs use `<method>/<dataset>/drop_out/drop_out_0.1/rep1.csv`, ratios 0.1–0.5, and repeats 1–3 corresponding to seeds 123–125. Means use finite genes within each run. The 15 scstGCN/Dataset 5 attempts were non-evaluable and are reported as such; they are not imputed. Paired changes use the same seed's 10% reference, 95% Student t intervals and two-sided paired 10%/50% endpoint tests. Holm adjustment is within each dataset–metric family across available methods. Intermediate ratios are descriptive. The raw-metric figures show means and sample SD with shared per-metric axes across methods within a dataset.

## 8. Usability and optional plotting interfaces

```bash
python -m pipeline.analysis.cli usability --input DATA/usability_observations.csv --output-dir RESULTS/usability
python -m pipeline.analysis.run_analysis --job configs/analysis_job.example.json
```

The usability input is a raw observation table with method names in its first column. It is not a table of already normalized scores. The JSON interface supports library functions for spatial expression examples and additional plotting. Mark path arguments with `{"__path__": "relative/path"}`, dataframe arguments with `{"__csv__": "relative.csv"}`, and images with `__image__`. Marked paths resolve relative to the job file. No private machine paths are needed.

## Metric file layout

The evaluators export CSV files into their requested output directories. Keep one CSV per method output, dataset, condition and computational repeat. Before `insample`, place or copy each baseline CSV into the following layout beneath `--metric-root`:

| Platform | Relative path |
| --- | --- |
| Xenium | `<method>/<dataset>/raw_resolution/<rep>/fixed_spotSize/spotDistance_100.csv` |
| Visium HD | `<method>/<dataset>/raw_resolution/<rep>/fixed_spotSize/spotDistance_96.csv` |
| CosMx | `<method>/<dataset>/<rep>.csv` |

Here `<rep>` is `rep1`, `rep2` or `rep3`. For example, copy the iSTAR Xenium baseline `rep1.csv` from the evaluator output into `DATA/metric_csv/iSTAR/Xe_hsa_breast_c/raw_resolution/rep1/fixed_spotSize/spotDistance_100.csv`. Keep the six metrics and gene identifiers unchanged. SuperST's two outputs use separate `SuperST_t` and `SuperST_t_g` directories. The evaluator's own output filename is not a simulation-condition identifier; assign the condition from the actual input used for that run. Resolution analyses use the analogous per-condition layout documented above.
