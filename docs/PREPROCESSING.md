# Preparing method inputs

Platform preprocessing remains separate from method execution. Use the retained
notebooks/scripts in `pipeline/preprocessing/` to register images and generate
the required simulated LR counts and matched fixed HR references.

Install `requirements-preprocessing.txt` and, for notebook execution,
`requirements-notebooks.txt`. Visium HD also supports CPU aggregation without
CuPy/Numba. R/HEARTSVG and interactive napari/CosMx registration use their own
environments.

## Platform preparation

- Xenium: align the H&E image using the exported registration matrix, filter and
  convert transcript coordinates, then generate pseudo-spot counts. Separate
  notebooks generate registered dropout inputs and FOV subsets.
- Visium HD: retain the actual study SVG list, rectify the image/coordinates
  using the supplied setup and aggregate 8-micrometre bins. The aggregation
  script optionally accepts `--gene_names` to restrict sparse gene rows before
  aggregation. Native H5 export is provided separately.
- CosMx: stitch/import exports, complete interactive registration, apply the
  actual test-region transform/crop and generate pseudo-spot counts. The
  registration notebook includes numerical coordinate-export functions. The retained
  workflow generates default circular pseudo-spots (55 µm diameter, 100 µm
  center spacing) and writes cnts.tsv and locs-raw.tsv in SlideOutputDir.
  CosMx is not included in the initial-resolution analysis.

Set the notebook's input paths, actual pixel sizes and finalized study bounds
before execution. The retained numeric examples are dataset-specific settings.
Confirm image/transcript alignment visually before method execution.

### CosMx registration and selected FOVs

`CosMx/1.Cosmx_napari_matrix_simulation.ipynb` retains the density-assisted
registration and numerical transcript-coordinate export. Use the
[manufacturer's alignment tutorial](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/napari-alignment-tools/)
for the interactive controls. Set `ALIGNMENT_FOV_RANGE` for the density reference
and `ANALYSIS_FOV_RANGE` for the transcripts to export. The examples `(1,36)`
and `(90,93)` demonstrate selecting subsets; they are not required ranges.
Different ranges must share the original global coordinate system and a
registration valid over the analysis region. Keep the density layer registered
to the CosMx reference when aligning H&E.

Set `COSMX_MPP` from the CosMx image metadata and `HE_MPP` from the supplied H&E
image. The example `HE_MPP=0.24891` is the pixel size of the study H&E, not a
platform constant. Cropping retains this pixel size; resizing requires an
updated value. The notebook saves `registration.json` and `pixel-size-raw.txt`.
`CosMx/3.different_formats.ipynb` reads the saved pixel size rather than using a
second hardcoded value. Crop bounds are explicit, and out-of-range transcripts
stop execution for review rather than being silently discarded.

The notebook exports `transcripts_use.tsv`, `he-raw.tiff`, `cnts.tsv` and
`locs-raw.tsv` for the selected region. Use the common packager below to assemble
complete model inputs. The auxiliary CosMx format notebook also exports padded
images; for a complete SuperST input, use the common SuperST preparation route
below instead of treating those image files alone as a complete input.

### Visium HD SVG input

To reproduce the study, use the finalized study SVG list. If generating a new
list with the retained HEARTSVG procedure, first export the native 8 µm bins:

```bash
python pipeline/preprocessing/VisiumHD/export_native_bins.py --matrix-dir DATA/DATASET/binned_outputs/square_008um --positions DATA/DATASET/binned_outputs/square_008um/spatial/tissue_positions.parquet --output-dir RESULTS/native_8um
```

Set `dataset_dir` in `1.selectSVG.ipynb` to this output and `result_dir` to the
desired gene-list directory. This exports native counts without coarsening or
reordering their barcodes; positions are matched by barcode. Writes are chunked,
but full-transcriptome bin-level TSVs can still be large. Do not run SVG selection
on coarsened simulated inputs. The resolution notebook requires the finalized
`GENE_NAMES_FILE` and passes it to aggregation without selecting genes again.
Raw clinical data and finalized per-dataset input lists are external data
resources, not inferred from prediction performance.

## Package simulated counts

The common input packager accepts an existing spot-by-gene matrix, corresponding
pixel coordinates and H&E image:

```bash
python -m pipeline.preprocessing.prepare_inputs --counts DATA/cnts.tsv --locations DATA/locs-raw.tsv --image DATA/he-raw.tiff --output-dir DATA/prepared --pixel-size 1.7 --spot-diameter 55 --dataset Xe_hsa_breast_c --gene-list DATA/gene-names.txt
```

Use the actual physical pixel size, diameter/side length and study gene list.
Pass `--radius-pixels` when the recorded method input uses a rounded radius.
Optional `--features` supplies a headered TSV with `id/name` metadata.

The packager aligns feature metadata with count columns, checks barcode
identities and writes matrix-format files plus the H5, position tables, images
and scale factors required by Visium-format methods. Counts and image coordinates
are retained. Both generated PNGs use the supplied image dimensions and scale
factor 1. Use an empty destination directory.

Prepare separate method/condition/repeat directories. SuperST image padding and
XFuse conversion/cropping use the explicit study settings described in
[method commands](METHODS.md). See [data requirements](DATA_REQUIREMENTS.md) for
input schemas and matched-ground-truth requirements.

## Connect platform outputs to method inputs

Each simulation condition produces a separate pair of counts and locations.
Use both from the same condition and the matching aligned H&E:

| Platform | Counts/locations | H&E |
| --- | --- | --- |
| Xenium baseline | `raw_resolution/matrix/fixed_spotSize/spotDistance_100/{cnts.tsv,locs-raw.tsv}` | Aligned dataset `he-raw.tiff` |
| Visium HD | The chosen `raw_resolution/matrix/<scheme>/<condition>/{cnts.tsv,locs-raw.tsv}` | `data_rectified/he_rectified.tiff`, matched to the adjusted positions |
| CosMx baseline | `SlideOutputDir/{cnts.tsv,locs-raw.tsv}` | `SlideOutputDir/he-raw.tiff` |

The Xenium baseline-format, FOV and dropout notebooks read the counts/locations
from the dataset root. Copy the two baseline files from the first row there
before those notebooks; keep the original aligned H&E and full HR transcript
reference. Do not copy another resolution condition into that baseline slot.
For each condition, pass the table's matching files to `prepare_inputs`. Native
H5-only notebook exports are alternatives for matrix conversion, not substitutes
for the images, scale metadata and gene list of a complete method input.

For SuperST, first run `prepare_inputs`, then set `prepared_dir` and a new
`superst_dir` in `methods/prepare_superst_image_size.ipynb`. It copies the complete
input and replaces only the spatial images with bottom/right-padded images.
Use the resulting directory with the launcher and the value in
`superst-target-size.txt`. Padding does not change the original coordinate origin,
counts or fixed HR reference. All other methods use their appropriate prepared
inputs as described in `METHODS.md`.

## Running a notebook with parameters

```bash
python -m pipeline.notebook_job --notebook pipeline/preprocessing/Xenium/2.Xenium_spot_matrix_simulate.ipynb --parameters configs/xenium_paths.example.json --output RESULTS/xenium.executed.ipynb --cwd .
```

Parameters replace existing top-level setup assignments. Unknown parameter names
are rejected; calculations inside functions are retained. The executed notebook
is written separately from the source. Complete interactive registration before
running dependent noninteractive preparation steps.

The Xenium example uses the package root as its kernel working directory, so
`DATA/` and `RESULTS/` in its JSON resolve there. The Visium HD matrix notebook
uses its own folder to locate neighboring helper scripts; open it from
`pipeline/preprocessing/VisiumHD/` and use absolute dataset paths, or set its
`SCRIPT_DIR` explicitly. Do not apply one working-directory rule to every
platform notebook. The CosMx registration notebook is interactive and must not
be executed unattended as a complete batch job.
