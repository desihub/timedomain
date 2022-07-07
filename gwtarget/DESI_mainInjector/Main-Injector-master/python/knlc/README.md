# KNLC

A tool to predict the brightness potential GW counterparts.

## Usage

As an example,

`python kn_brightness_estimate.py --distance 100 --distance_err 20 --time_delay 6 --magplot_file kn_magplot.png --expplot_file kn_expplot.png --report_file kn_report`

will produce brightness estimates for a merger at 100 +/- 20 Mpc at a time 6 hours post-merger.
Plots for the Main-Injector website will be saved in the current working directory in `kn_magplot.png` and `kn_expplot.png`.
A human-readable report will be saved to `kn_report.txt`.
This report will contain the expected magnitude for GW170817.
A csv file of the mags in each band at each percentile threshold will be saved to `kn_report.csv`

At present, only `time_delay` values of less than 56.3 hours are implemented, but this is an active area of developement.
