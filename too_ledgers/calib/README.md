Ledger construction using data from DECAT and the SMBBH list in DESI
calibration tiles. These are targets for DESI calibration fields, located here:

https://desi.lbl.gov/trac/wiki/SurveyOps/CalibrationFields

DECAT data were selected from Rob Knop's pipeline at https://decat-webap.lbl.gov/decatview.py/ using the following settings:
* Use real/bogus type: 2
* Use Version tag: latest (only use this; version tags aren't implemented yet)
* Search only proposal IDs: select just 2019A-0065 (but note: that's only for the initial search!)
* Search objects detected starting: box checked, value 2022-10-01
* Search objects detected ending: box checked, value 2022-10-31
* S/N cut: box not checked
* Only high r/b objects in initial search: box checked
* Filter detection counts starting: box checked, value 2022-10-01
* Filter detection counts ending: 2022-10-31
Uncheck everything down to:
* Seen in at least this many bands: box checked, value 3
* Min # detections: box unchecked
* Min # high r/b detections: box checked, value 6
* Max # detections outside date range: box unchecked
* Max high r/b detections outside date range: box checked, value 2
Can leave latitude limits unchecked.
