# From Lei, Jan 20, 2023

# ****************************** Definition of columns in the table DECam_COSMOS_*.csv ****************************** #
#
# [1] objid: Object Unique ID, e.g., C202210060223125m045817
#            objid with starting 'A': the separation to nearest (Legacy Survey DR9) extended source < 0.3"
#            objid with starting 'C': the separation to nearest (Legacy Survey DR9) extended source between 0.3" - 1.0"
#            objid with starting 'T': the separation to nearest (Legacy Survey DR9) extended source > 1.0"
#            P.S. the objects with starting 'A' are possibly AGNs.
#
# [2] skycoord_obj: Object Sky Coordinate (string), e.g., 02:23:12.565 -04:58:17.903
#     * Remarks: it is the coordinate of the first alert of this object.        
#                    
# [3] ra_obj: Object Sky Coordinate (RA, deg), e.g., 35.80235817
#
# [4] dec_obj: Object Sky Coordinate (DEC, deg), e.g., -4.971639968
#
# [5] ra_obj_maxSNR: Refined Object Sky Coordinate (RA, deg), e.g., 35.80239381
#     * Remarks: it is the coordinate of the alert with maximum SExtractor SNR_WIN of this object.
#                the force photometry is based on such refined coordinate.
#
# [6] dec_obj_maxSNR: Refined Object Sky Coordinate (DEC, deg), e.g., -4.971643086
#
# [7] num_alert: The total number of alerts, e.g., 10
#     * Remarks: a SExtractor detection is an alert iff it passes the CNN stamp classifier 
#                and its force photometry has a SNR > 5.0
# 
# [8] date_first_alert: The UTC time of the first alert, e.g., 2022-10-06T05:32:39.469
#
# [9] mag_first_alert: The magnitude of the first alert, e.g., 22.05587399
#
# [10] filter_first_alert: The filter (g/r/i for XMM fields) of the first alert, e.g., i
#
# [11] date_peak_alert: The UTC time of the alert with peak magnitude, e.g., 2022-10-06T05:32:39.469
#
# [12] mag_peak_alert: The magnitude of the alert with peak magnitude, e.g., 22.05587399
#
# [13] filter_peak_alert: The filter (g/r/i for XMM fields) of the alert with peak magnitude, e.g., i
#
# [14] date_last_alert: The UTC time of the last (until Nov 24) alert, e.g., 2022-11-25T05:00:25.433
#
# [15] mag_last_alert: The magnitude of the last (until Nov 24) alert, e.g., 23.07915265
#
# [16] filter_last_alert: The filter (g/r/i for XMM fields) of the last (until Nov 24) alert, e.g., r
#
# [17] lsdr9_id.pseudo_host: pseufo host galaxy means the nearest extended source in Legacy Survey DR9, e.g., 0357m050_6582
#      * Remarks: lsdr9_id has a format like legacy_survey_dr9_tractor_name + '_' + objid_in_tractor_catalog
# 
# [18] lsdr9_ra.pseudo_host: pseufo host galaxy ra in Legacy Survey DR9, e.g., 35.80260165
#
# [19] lsdr9_dec.pseudo_host: pseufo host galaxy dec in Legacy Survey DR9, e.g., -4.971537074
#
# [20] lsdr9_gmag.pseudo_host: pseufo host galaxy g-band magnitude in Legacy Survey DR9, e.g., 19.1067341
#
# [21] lsdr9_rmag.pseudo_host: pseufo host galaxy r-band magnitude in Legacy Survey DR9, e.g., 18.07562858
#
# [22] lsdr9_zmag.pseudo_host: pseufo host galaxy z-band magnitude in Legacy Survey DR9, e.g., 17.29460895
#
# [23] z.pseudo_host: redshift of the pseufo host galaxy, e.g., 0.149849996
#      * Remarks: searched in GAMA-DR3, GLADE, SDSS, SixdF, GWGC, NED within 1".
#
# [24] e_z.pseudo_host: redshift error of the pseufo host galaxy, e.g., NaN
#
# [25] z_type.pseudo_host: redshift type, spectroscopic or photometric, e.g., spec
#
# [26] z_source.pseudo_host: redshift source, e.g., GAMA-DR3
#
# NOTE: the finding chart of each candidate has a filename format such as 
#       10_C202210060223125m045817.pdf, which is [num_alerts]_[objid].pdf
#