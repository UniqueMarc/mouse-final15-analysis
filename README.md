# Mouse bioluminescence image analysis

Quantify three mouse regions per image across 15 TIFF images: five treatment groups at 2, 4 and 8 hours, giving 45 ROI measurements.

## Install and run

Use Python 3.12. From the project folder:

```sh
python -m pip install -r requirements.txt
python analyze_mice.py
```

To choose an output folder and save ROI overlays:

```sh
python analyze_mice.py --output-dir analysis_output --write-qc
```

## Inputs

- `data/images/`: 15 input TIFF images.
- `data/manifest.json`: image paths, SHA-256 checksums, color-bar limits and ROI coordinates.
- `data/image_index.csv`: image names, time points and treatment groups.

Keep the folder structure intact. ROI coordinates use exclusive upper bounds. The program verifies image checksums before quantification.

## Results

Results are written to `outputs/` by default:

| File | Contents |
| --- | --- |
| `quantification.csv` | One row per ROI, including the image-derived integrated signal |
| `summary.csv` | Mean and sample standard deviation of three ROIs per image |
| `verification.json` | Input checksums, measurement definitions and output counts |
| `roi_qc/` | ROI overlays when `--write-qc` is enabled |

Open CSV files in Excel or another analysis program. Supplied result tables are in `results/`.

## Interpret the measurements

The program matches colored ROI pixels to the color bar in the same image and sums the mapped values. Results are **image-derived integrated bioluminescence signal (a.u.)**. No background subtraction, ROI-area normalization or cross-image normalization is applied. A zero means no pixels passed the color-selection criteria.

The inputs are exported RGB images; this signal is not an instrument-calibrated photon flux. ROI positions identify locations within an image and do not establish animal identities across time points. The program reports descriptive summaries without group comparisons or repeated-measures tests.

## Check the installation

```sh
python -m unittest discover -s tests -v
```
