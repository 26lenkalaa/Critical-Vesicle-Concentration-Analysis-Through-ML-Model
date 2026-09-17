# CVC Absorbance Spectrum Analysis

This project determines Critical Vesicle/Micelle Concentration (CVC) from plate-reader
absorbance spectra using a ratiometric dye whose high/low wavelength absorbance ratio
shifts once vesicles or micelles start forming. What follows is not just a description
of the current pipeline but a record of how it got here, since most of what changed
along the way was driven by discovering that an earlier assumption did not match the
actual data, and that history is worth keeping rather than quietly overwriting.

## Where this started

The original two scripts were written during a summer of research in Dr. Maurer's lab
the summer after sophomore year of high school, as a first attempt at replacing manual
peak-picking with something automatic. Returning to the project years later to improve
it meant confronting assumptions that had gone untested since then, and most of what
follows is the record of finding and correcting those, not a rewrite from scratch.

The original two scripts, `sample15_absorbance_analysis.py` and
`sample26_absorbance_analysis.py`, worked like this: for each sample column, find that
sample's own peak absorbance independently in a high-wavelength window and a
low-wavelength window, baseline-subtract each against the absorbance at 600 nm, and take
the ratio. The last three columns in the file were assumed to be standards at known
concentrations of 0, 5, and 10, used to fit a linear calibration curve via
`numpy.polyfit`, which was then used to back-calculate a concentration for every other
column.

This did not work. Testing against the real data showed the last three columns in
Sample15 all returned a ratio of exactly 1.0, since none of them showed a distinct peak,
which flattened the calibration line to a slope of zero and made every back-calculated
concentration meaningless. Sample26's last three columns returned non-monotonic ratios,
producing a negative slope and negative estimated concentrations for real samples. That
was the first sign that "the last three columns are the standards" was a placeholder
assumption, not a confirmed fact about the data.

## First round of changes: automation without knowing the ground truth

The first rewrite (`automatedAnalysis.py`) tried to remove the standard-column
assumption by treating each file as one full serial dilution series and applying two
independent unsupervised breakpoint-detection methods: K-means clustering on the ratio
values, and a four-parameter sigmoid curve fit of ratio against log-concentration. Both
were printed side by side so they could be sanity-checked against each other. This was a
reasonable general-purpose approach given no other information, but it still rested on
assumptions that had not been confirmed: a single dilution factor, a single top
concentration, and every column belonging to one continuous series.

## The plate map: finding out each file wasn't one series at all

A photo of the lab's own well-assignment notes revealed the actual structure: each CSV
file only captured two rows of a much larger 96-well plate, and those two rows crossed
the boundaries of several differently-named sample groups, not one continuous dilution
series. Sample26 (rows A-B) actually contained three groups: seven replicate wells of a
single condition, one clean 10-well dilution series, and a partial, incomplete series cut
off by the edge of the uploaded rows. Sample15 (rows G-H) contained a similar mix, plus a
group that stayed flat across all nine of its wells with no transition at all, which
turned out to be a baseline/control group by design rather than a failure of the method.

This reframed the entire project. Reconstructing the plate map and re-running the
analysis by group, instead of by file, was what finally produced a clean, monotonic
dilution series where the earlier "last three columns" approach had produced nonsense.
It also surfaced a real anomaly: two wells in the middle of an otherwise low, flat series
read like high-concentration wells, which persisted under multiple ratio-calculation
methods and so was flagged as a genuine plate issue rather than a computational one.

## Fixing the ratio calculation itself

The original per-well method, where each column searched its own high and low windows
independently for a peak, turned out to have a quiet failure mode: for any well with no
real spectral shift, that search had nothing real to find, and noise pushed the "peak"
to wherever it happened to land, which collapsed the ratio to exactly 1.0000 for nearly
every below-CVC well. That looked like a real, clean baseline, but it was an artifact.
Switching to one fixed, shared high and low wavelength pair, chosen from where real peaks
consistently land in wells that do show a shift, replaced that artificial floor with a
proper continuous signal and made the transition curve visibly cleaner.

## Matching the actual breakpoint method used in the lab

A second whiteboard photo showed the deviation-from-linearity method actually used to
call CVC for this project: fit a line to the lowest-concentration points, extend it one
point at a time, and track R-squared. While R-squared stays high, those points are still
in the pre-transition linear region; the first point that causes a meaningful drop marks
the transition. This replaced the K-means/sigmoid approach as the primary method, since
those were reasonable unsupervised techniques but not what the lab's own methodology
actually was. One useful property of this method is that R-squared for a linear fit does
not change under a rescaling of evenly-spaced x-values, which meant the breakpoint could
be correctly located by well position even before the real concentration values (top
standard concentration, dilution spacing) were ever confirmed.

## Generalizing into a reusable pipeline

The project was then split into a reusable core (`cvc_core.py`) and a command-line
entry point (`analyze.py`), so that a new file could be dropped in without hand-editing
analysis code for it. This added a few things the earlier one-off scripts did not have:
automatic wavelength selection (picking whichever wavelength shows the most spread across
wells, rather than a hardcoded value), a declarative YAML plate-map format so a known
layout can be described in a few lines, group typing (`series`, `replicate`, or
`baseline`) so replicate and control wells stop being forced through breakpoint fitting,
and a statistical auto-segmentation fallback for files with no known plate map at all.

## Automating the pipeline and adding an AI review step

The pipeline was then wired into GitHub Actions, so pushing a new spectrum CSV into
`data/` triggers the full pipeline automatically and commits the results back to the
repository, with no manual steps required. An optional layer was added that sends only
derived summary statistics, never raw spectra, to the Claude API to flag anything
statistically suspicious and give a plain-language read of each result. This step is
deliberately non-blocking: if no API key is configured, the deterministic analysis still
runs and completes normally.

## The supervised-learning experiment

Most recently, a supervised classifier was added as an explicit experiment rather than a
replacement for the R-squared method. Using the two confirmed dilution series as labeled
training data (each well tagged above or below its confirmed breakpoint), a logistic
regression was trained to predict that label from the ratio and well position. Given
only twenty labeled wells total, leave-one-out cross-validation was used instead of a
standard train/test split, and came out around 80% accurate. Applying the trained model
back onto its own two training series, as a sanity check rather than real validation,
showed it did not reproduce either series' confirmed breakpoint exactly, which is the
expected and honest result of pooling two series of different noise levels into one
small model rather than a sign the code is broken. The path to a better version of this
is labeling more confirmed series as they become available, not tuning harder on these
twenty points.

## Where things stand now

Two dilution series across the two uploaded files are complete and trustworthy:
`original_concentrations_highprob_replace_smi5` (Sample26), with a high-confidence
breakpoint, and `replace_smi5_fm_finetune_concentrations` (Sample15), with a
lower-confidence breakpoint due to noise at its low-concentration end. Everything else in
the two files is either replicate-only, baseline-only by design, or an incomplete series
missing wells from rows that were never uploaded.

Two things remain unresolved by design rather than by oversight: the real top standard
concentration and dilution spacing for either plate have never been confirmed, so every
result is reported as "between well X and well Y" rather than a concentration in mM, and
the auto-segmentation fallback is a heuristic for when no plate map exists, not a
substitute for one.

## Requirements

```
pip install -r requirements.txt
```

Core dependencies: pandas, numpy, matplotlib, scikit-learn, scipy, pyyaml. The optional
AI review step additionally requires the `anthropic` package and an `ANTHROPIC_API_KEY`.

## Usage

```
python analyze.py data/your_file.csv --plate-map plate_maps/cvc_project.yaml --use-ai
python analyze.py data/your_file.csv --use-ai   # no known plate map, auto-segmented
```

## Files

- `analyze.py`, `cvc_core.py` — the current pipeline: cleaning, ratio calculation,
  grouping, R-squared breakpoint detection, optional AI review
- `plate_maps/cvc_project.yaml` — the confirmed plate layout for this project's
  formulation-comparison plate
- `automatedAnalysis.py` (supervised-learning version) — the logistic regression
  experiment described above, kept separate from the main pipeline since it is a proof
  of concept, not a production method
- `legacy/` — every earlier version of this pipeline, kept for history: the first
  per-file scripts, the first draft scripts, and the original K-means/sigmoid automation
  attempt
- `data/` — input CSVs
- `results/` — generated reports and plots