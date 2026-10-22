# Histological FOV similarity tool

This standalone tool extracts UNI histological features and searches for image regions similar to reference training FOVs. Similarity is an image-feature comparison, not a validated guarantee of super-resolution prediction accuracy.

```bash
python -m pip install -r tools/fov_similarity/requirements.txt
python tools/fov_similarity/extract_single_region_feature.py --help
python tools/fov_similarity/extract_top3_uni_features.py --help
python tools/fov_similarity/search_similar_regions.py --help
```

- `extract_top3_uni_features.py`: build reference features from training FOVs selected by the supplied whole-slide prediction-performance records.
- `extract_single_region_feature.py`: extract features for specified image bounds.
- `search_similar_regions.py`: compare reference features against windows in a user H&E image and export ranked windows, features and preview images.

Provide the actual UNI checkpoint and one reference `.npy` feature vector per search. Both are external resources. Reference and query features must use the same checkpoint, preprocessing, coordinate convention and image pixel-size assumptions. The scripts' `--help` lists the exact required input and output arguments. Do not add model weights or authentication tokens to this repository.

## Prepare and use references

Run the FOV analysis once for each dataset/method and write its output to `RESULTS/fov/<dataset>/<method>/`. The resulting `fov_gene_metrics.csv` contains the coordinates and whole-slide performance used to select up to three valid reference training FOVs. Numeric `1500` and labels such as `fov1500` are accepted. Regions without finite mean PARnorm are not selected.

```bash
python tools/fov_similarity/extract_top3_uni_features.py --base_csv_dir RESULTS/fov --he_base_dir DATA --output_base_dir RESULTS/reference_features --ckpt_path WEIGHTS/uni.bin --fov_size_target 1500
```

Use `--csv_filename fov_all.csv` for an existing collection with that filename. Each exported `.npy` is a single reference vector of shape `(D,)`. The search script accepts `(D,)` or `(1,D)`, not a stacked `(N,D)` bank. Run a separate search for each reference of interest; reference vectors are not averaged.
