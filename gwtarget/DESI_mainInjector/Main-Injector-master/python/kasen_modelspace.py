import os
import numpy as np
from knlc import kn_brightness_estimate

def run_ap_mag_for_kasen_models (filter, distance, dist_err, days_since_burst, 
        kasen_fraction, data_dir="./", fast=False, doPlots=True) :
    report_file = os.path.join(data_dir, "kn_report")
    if not fast :
        knlc_dir = os.path.join(os.getenv("DESGW_DIR", data_dir), "knlc")
        print(knlc_dir)
        code = os.path.join(knlc_dir, "kn_brightness_estimate.py")
        cmd = "python {} --distance {} --distance_err {} --time_delay {} ".format(
            code, distance, dist_err, days_since_burst)
        cmd = cmd + "--fraction {} ".format(kasen_fraction)
        if doPlots :
            cmd = cmd + "--magplot_file kn_mag_plot.png "
            cmd = cmd + "--expplot_file kn_exp_plot.png "
            cmd = cmd + "--report_file {} ".format(report_file)
        os.system(cmd)

        file = report_file + ".txt"
        fd = open(file,"r")
        for i in range(0,16): #fd.readline()
            line = fd.readline().split()
            apparent_mag = dict()
            try:
                apparent_mag["g"] = np.float(line[0])
                apparent_mag["r"] = np.float(line[3])
                apparent_mag["i"] = np.float(line[6])
                apparent_mag["z"] = np.float(line[9])
                ap_mag = apparent_mag[filter]
            except:
                pass

    else :
        kn_calc         = kn_brightness_estimate.KNCalc(distance, dist_err, days_since_burst)
        percentile_dict = kn_brightness_estimate.calc_mag_fractions(kn_calc.template_df_full)
        cutoffs         = kn_brightness_estimate.mags_of_percentile(
            kasen_fraction, percentile_dict)
        kn_brightness_estimate.make_output_csv(
            np.linspace(0., 100., 101), 
            percentile_dict, 
            write_answer=True, 
            flt=filter, 
            fraction=kasen_fraction,
            data_dir=data_dir)

        file = os.path.join(data_dir, "answer_{}.txt".format(filter))
        ap_mag = np.genfromtxt(file)
 
    return ap_mag
