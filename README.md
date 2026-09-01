# Absorbance Spectrum Analysis

This project analyzes plate-reader absorbance spectra exported as CSV files (`Sample15_Absorbance_Spectrum.csv`, `Sample26_Absorbance_Spectrum.csv`) and estimates sample concentrations from a standard curve. Each sample column in a file is treated independently: the script finds that sample's own absorbance peak in a high-wavelength region and its own peak in a low-wavelength region, baseline-subtracts each against the absorbance at 600 nm, and takes the ratio of the two. A handful of columns with known concentrations are then used to fit a straight-line calibration curve of concentration versus ratio, and that line is used to back-calculate an estimated concentration for every other sample on the plate.

There are two scripts, `sample15_absorbance_analysis.py` and `sample26_absorbance_analysis.py`, kept deliberately separate rather than merged into one shared module, since the two plates are handled as independent runs with slightly different entry points (Sample15 prompts for a filename and stock concentration at the command line; Sample26 runs directly against its known filename).

## Input file format

The CSVs exported by the plate reader are not plain data tables. Each file starts with several lines of instrument metadata (software name, run date, wavelength range, plate number) before the actual header row, which begins with `Wavelength,`. The scripts locate that header row automatically by scanning for a line that starts with `Wavelength,` specifically, rather than assuming a fixed number of rows to skip, since an earlier line in the file (`Wavelength: 400-600 nm`) also happens to start with the word "Wavelength" and will produce a parse error if matched carelessly.

Below the header row, each row is one wavelength (typically 400 to 600 nm, one row per nanometer) and each column after `Wavelength` is one sample well, labeled with an ID like `Un0001 (G01)` or `Blank1 (H10)`. The scripts read this shape directly from the file rather than hardcoding row or column counts, so a plate with a different wavelength range or a different number of wells will still load correctly.

## Methodology

For each sample column, the script searches two wavelength windows independently: 400–600 nm for the "high peak" and 400–530 nm for the "low peak." Each peak value is baseline-subtracted using that same sample's absorbance at 600 nm, and the ratio of the high-peak subtraction to the low-peak subtraction is recorded as that sample's signal.

A set of columns with known concentrations (by default, the last three columns in the file, corresponding to concentrations of 0, 5, and 10) are used to fit a linear calibration curve of concentration against this ratio via `numpy.polyfit`. The resulting slope and intercept are then used to back-calculate an estimated concentration for every remaining sample column.

## Known limitation: standard column selection is not yet confirmed

Testing against the actual uploaded data showed that treating "the last three columns in the file" as the known standards does not reliably work for either plate. In Sample15, the last three columns (`Un0020`, `Un0021`, `Blank1`) all return a ratio of exactly 1.0, since none of them show a distinct peak in the 550–600 nm region, which flattens the calibration line to a slope of zero and makes every back-calculated concentration meaningless. In Sample26, the last three columns return non-monotonic ratios (2.51, 1.0, 1.0), producing a negative slope and negative estimated concentrations for most real samples. Sample26's column headers are also not in strict numeric order in the raw file (`Un0019`–`Un0021` appear before `Un0010`–`Un0018`), so column position alone is not a safe way to identify the standards.

Until the actual standard well IDs are confirmed for each file, the `standard_concentrations` and `n_standards` arguments in `main()` should be treated as placeholders, and any concentration output from the current scripts should not be relied on for a report. The fix, once the correct wells are known, is to select standards by well label instead of column position.

## Requirements

- Python 3.9+
- pandas
- numpy

Install with:

```
pip install pandas numpy
```

## Usage

Sample15 prompts for the CSV filename and, optionally, stock solution details at the command line:

```
python sample15_absorbance_analysis.py
```

Sample26 runs directly against its expected filename in the working directory:

```
python sample26_absorbance_analysis.py
```

Each script prints, per sample column: the high peak value and wavelength, the low peak value and wavelength, and the high/low ratio. It then prints the fitted calibration curve (slope and intercept) and an estimated concentration for every non-standard sample.

## Files

- `sample15_absorbance_analysis.py` — analysis script for Sample15, with an interactive filename and stock-concentration prompt
- `sample26_absorbance_analysis.py` — analysis script for Sample26, run against its known filename
- `Sample15_Absorbance_Spectrum.csv`, `Sample26_Absorbance_Spectrum.csv` — raw plate reader exports

## Next steps

- Confirm the actual standard well IDs and their known concentrations for each file, and update the scripts to select standards by well label rather than column position
- Decide whether `get_concentration_info()` in the Sample15 script (which computes a stock dilution concentration) should feed into the calibration curve's known concentrations, or remain a separate, informational calculation
- Once the standards are fixed, sanity-check that the calibration curve is monotonic (absorbance ratio should increase or decrease consistently with concentration) before trusting any back-calculated sample concentration