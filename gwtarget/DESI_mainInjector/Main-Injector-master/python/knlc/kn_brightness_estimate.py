# A module to print brightness estimates

import matplotlib
matplotlib.use("Agg");
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from optparse import OptionParser
import pandas as pd
from scipy.stats import norm
import sys
import os

from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value
import astropy.units as u

# Handle command-line arguments
class KNCalc():
    def __init__(self, distance, distance_err, time_delay):
        
        self.distance = distance
        self.distance_err = distance_err
        
        # Convert time_delay from hours to days
        if float(time_delay) > 400.8:
            print("Currently, only time delays less than 400.8 hours (16.7 days) post merger are supported")
            sys.exit()
        self.delta_mjd = round(float(time_delay) / 24.0, 1)
        
        # Set directory for lookup table
        knlc_dir = os.getenv("DESGW_DIR", "./")
        if knlc_dir != "./" : knlc_dir = knlc_dir + "/knlc/"
            
        # Choose lookup table based on time_delay
        if self.delta_mjd < 2.3:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry.csv')
        elif self.delta_mjd < 4.7:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_2.csv')
        elif self.delta_mjd < 7.1:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_3.csv')
        elif self.delta_mjd < 9.5:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_4.csv')
        elif self.delta_mjd < 11.9:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_5.csv')
        elif self.delta_mjd < 14.3:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_6.csv')
        elif self.delta_mjd < 16.7:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_7.csv')
        df['ZMEAN'] = np.mean(df[['ZMIN', 'ZMAX']].values, axis=1)

        # Mean distance calculation 
        mean_z = z_at_value(cosmo.luminosity_distance, float(self.distance) * u.Mpc)
        template_df_mean = df[(df['ZMIN'].values < mean_z) & (df['ZMAX'].values > mean_z) & (df['DELTA_MJD'].values == self.delta_mjd)].copy().reset_index(drop=True)
        template_df_mean['WEIGHT'] = 1.0 / template_df_mean.shape[0]
        self.template_df_mean = template_df_mean

        # Full distance calculation
        template_df_full = df[df['DELTA_MJD'].values == self.delta_mjd].copy().reset_index(drop=True)
        weights = [norm.pdf(x.value, loc=float(self.distance), scale=float(self.distance_err)) for x in cosmo.luminosity_distance(template_df_full['ZMEAN'].values)]
        template_df_full['WEIGHT'] = weights / np.sum(weights)
        self.template_df_full = template_df_full

        return

        

### Functions to calculate metrics of interest
def weighted_average(quantity, weights):
    return np.dot(quantity, weights) / np.sum(weights)

def get_all_mags(data):
    out_dict = {'g_mag': weighted_average(data['MAG_g'].values, data['WEIGHT'].values),
                'r_mag': weighted_average(data['MAG_r'].values, data['WEIGHT'].values),
                'i_mag': weighted_average(data['MAG_i'].values, data['WEIGHT'].values),
                'z_mag': weighted_average(data['MAG_z'].values, data['WEIGHT'].values),
                'g_magerr': weighted_average(data['MAGERR_g'].values, data['WEIGHT'].values),
                'r_magerr': weighted_average(data['MAGERR_r'].values, data['WEIGHT'].values),
                'i_magerr': weighted_average(data['MAGERR_i'].values, data['WEIGHT'].values),
                'z_magerr': weighted_average(data['MAGERR_z'].values, data['WEIGHT'].values)}
    return out_dict

def gw170817(data):
    blue_df = data[data['SIM_TEMPLATE_INDEX'].values == 178.0]
    red_df = data[data['SIM_TEMPLATE_INDEX'].values == 224.0]

    blue_mags = get_all_mags(blue_df)
    red_mags = get_all_mags(red_df)

    return blue_mags, red_mags

def print_dict(d, title='', outfile=None):

    outdata = [title + '\n--------------------------------------------------------------\n',
               "\tg\t\tr\t\ti\t\tz\n",
               "--------------\t--------------\t--------------\t--------------\t\n", 
               "%.2f +/- %.2f\t%.2f +/- %.2f\t%.2f +/- %.2f\t%.2f +/- %.2f\n\n" %(d['g_mag'], d['g_magerr'],
                                                                                d['r_mag'], d['r_magerr'],
                                                                                d['i_mag'], d['i_magerr'],
                                                                                d['z_mag'], d['z_magerr'])]
    if outfile:
        stream = open(outfile + '.txt', 'a+')
        stream.writelines(outdata)
        stream.close()
    else:
        for row in outdata:
            print((row[:-1]))

    return


def calc_mag_fractions(data):

    averaged_data = pd.DataFrame(data=[get_all_mags(data[data['SIM_TEMPLATE_INDEX'].values == y]) for y in np.unique(data['SIM_TEMPLATE_INDEX'].values)])
    averaged_data['SIM_TEMPLATE_INDEX'] = np.unique(data['SIM_TEMPLATE_INDEX'].values)

    percentile_levels = np.linspace(0.0, 100.0, 101)
    percentile_dict = {}
    for band in ['g', 'r', 'i', 'z']:
        percentile_dict['%s_cutoff' %band] = np.nanpercentile(averaged_data['%s_mag' %band].values, q=percentile_levels)

    return percentile_dict


def mags_of_percentile(cutoff, percentile_dict):
    if cutoff < 1.0:
        cutoff *= 100

    index = int(round(cutoff))
    return {band: percentile_dict['%s_cutoff' %band][index] for band in ['g', 'r', 'i', 'z']}

def make_output_csv(cutoffs, percentile_dict, outfile=None, return_df=False, write_answer=False, flt='', fraction=90.0):

    out_data = [mags_of_percentile(cutoff, percentile_dict) for cutoff in cutoffs]
    out_df = pd.DataFrame(out_data)
    out_df['PERCENTILE'] = cutoffs
    if outfile:
        out_df.to_csv(outfile + '.csv', index=False)

    if write_answer:
        if fraction < 1.0:
            farction *= 100
        #get closest index to fraction
        closest_index = np.argmin(np.abs(float(fraction) - out_df['PERCENTILE'].values))
        stream = open('answer_%s.txt' %flt, 'w+')
        stream.write('%.2f' %out_df[flt].values[closest_index])
        stream.close()
    
    if return_df:
        return out_df


def make_plot(percentile_dict, blue, red, title='', outfile=None, fraction=None):
    plt.figure()

    percentile_levels = np.arange(len(percentile_dict['g_cutoff']))

    color_dict = {'i': 'brown', 'g': 'green', 'r': 'red', 'z': 'dimgray'}

    interps = {}
    for band in ['g', 'r', 'i', 'z']:
        m0=get_m0(band)
        plt.plot(percentile_dict['%s_cutoff' %band], percentile_levels, lw=2, label=band, color=color_dict[band])
        plt.axvline(x=m0,color=color_dict[band], lw=1, ls=':')
        interps[band] = interp1d(percentile_dict['%s_cutoff' %band], percentile_levels)
        plt.errorbar(blue['%s_mag' %band], interps[band](blue['%s_mag' %band]), xerr=blue['%s_magerr' %band], 
                     capsize=2, marker='o', markerfacecolor='None', markeredgecolor=color_dict[band], color=color_dict[band])
        plt.errorbar(red['%s_mag' %band], interps[band](red['%s_mag' %band]), xerr=red['%s_magerr' %band], 
                     capsize=2, marker='s', markerfacecolor=color_dict[band], markeredgecolor=color_dict[band], color=color_dict[band])
    
        if band == 'z':
            #plt.errorbar([22.5], [25], xerr=[0.3], capsize=2, marker='o', markerfacecolor='None', markeredgecolor=color_dict[band], color=color_dict[band])
            #plt.text(23, 25, 'GW170817-blue', verticalalignment='center', fontsize=12)
            #plt.errorbar([22.5], [15], xerr=[0.3], capsize=2, marker='s', markerfacecolor=color_dict[band], markeredgecolor=color_dict[band], color=color_dict[band])
            #plt.text(23, 15, 'GW170817-red', verticalalignment='center', fontsize=12)
            plt.errorbar([18.6], [60], xerr=[0.3], capsize=2, marker='o', markerfacecolor='None', markeredgecolor='k', color='k')
            plt.text(19.1, 60, 'GW170817-blue', verticalalignment='center', fontsize=12)
            plt.errorbar([18.6], [50], xerr=[0.3], capsize=2, marker='s', markerfacecolor='k', markeredgecolor='k', color='k')
            plt.text(19.1, 50, 'GW170817-red', verticalalignment='center', fontsize=12)
            plt.plot([18.3,18.9], [40,40], lw=1, ls=':', color='k')
            plt.text(19.1, 40, '90s 10$\sigma$ limit', verticalalignment='center', fontsize=12)
    
            plt.xlabel("magnitude", fontsize=14)
            plt.ylabel("percent (<magnitude) ", fontsize=14)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)

    if fraction:
        fraction = float(fraction)
        if fraction < 1.0:
            fraction *= 100
        #plt.axhline(y=fraction, color='black', lw=1, ls='--')

    #plt.grid()
    plt.legend(fontsize=14, loc='upper left', frameon=False)
    plt.ylim(0,100)
    plt.xlim(18,26)
    plt.title(title, fontsize=14)
    
    if outfile:
        plt.savefig(outfile)

    plt.close()

    return

def get_percentile_at_exp_time(exptime, band, percentile_dict):
    m0 = get_m0(band)
    mag = m0 + 1.25 * np.log10(exptime / 90.0)

    df = make_output_csv(np.linspace(0.0, 100.0, 101), percentile_dict, outfile=None, return_df=True)

    index_of_closest_mag = np.argmin(np.abs(mag - df[band].values))
    percentile = df.iloc[index_of_closest_mag]['PERCENTILE']
    return percentile


def get_exptime(m0, mag):
    return 90.0 * 10 ** ((mag - m0) / 1.25)

def get_m0(band):
    out_dict = {'g': 23.4, 'r': 23.1, 'i': 22.5, 'z': 21.8, 'Y': 20.3}
    return out_dict[band]

def make_exptime_plot(percentile_dict, title='', outfile=None):
    percentile_levels = np.arange(len(percentile_dict['g_cutoff']))

    color_dict = {'i': 'brown', 'g': 'green', 'r': 'red', 'z': 'dimgray'}

    #color_dict = {'g': 'darkblue', 'r': 'darkgreen', 'i': 'darkred', 'z': 'black'}

    fig,ax1 = plt.subplots()
#    ax2 = ax1.twiny()
#    ax3 = ax1.twinx()

    for band in ['g', 'r', 'i', 'z']:
        m0 = get_m0(band)
        exptimes = get_exptime(m0, percentile_dict['%s_cutoff' %band])
        #plt.plot(exptimes, percentile_levels, lw=2, label=band, color=color_dict[band])
        ax1.plot(exptimes, percentile_levels, lw=2, label=band, color=color_dict[band])
#        ax2.plot(percentile_dict['%s_cutoff' %band], percentile_levels, lw=2, ls='--', color=color_dict[band])
#        ax3.plot(percentile_dict['%s_cutoff' %band], percentile_levels, lw=2, ls='--', color=color_dict[band])

    ax1.axvline(x=90.0, ls='--', lw=0.5)

    percentiles_at_90_sec = {band: get_percentile_at_exp_time(90.0, band, percentile_dict) for band in ['g', 'r', 'i', 'z']}
    

    ax1.text(105, 25, "90 sec in g = %.1f %%" %(percentiles_at_90_sec['g']), horizontalalignment='left', verticalalignment='center', fontsize=11)
    ax1.text(105, 20, "90 sec in r = %.1f %%" %(percentiles_at_90_sec['r']), horizontalalignment='left', verticalalignment='center', fontsize=11)
    ax1.text(105, 15, "90 sec in i = %.1f %%" %(percentiles_at_90_sec['i']), horizontalalignment='left', verticalalignment='center', fontsize=11)
    ax1.text(105, 10, "90 sec in z = %.1f %%" %(percentiles_at_90_sec['z']), horizontalalignment='left', verticalalignment='center', fontsize=11)

    ax1.set_xlabel("exptime * teff", fontsize=14)
    ax1.set_ylabel("percent (< 10$\sigma$ limiting mag)", fontsize=14)
#    ax3.set_ylabel("percent (< mag)", fontsize=14)
#    ax2.set_xlabel("mag", fontsize=14)
    #ax1.xticks(fontsize=12)
    #ax1.yticks(fontsize=12)
    ax1.set_xlim(30,300)
    ax1.set_ylim(0,100)

    #plt.xlabel("exptime * teff", fontsize=14)
    #plt.ylabel("percent (<maglim10)", fontsize=14)
    #plt.xticks(fontsize=12)
    #plt.yticks(fontsize=12)
    #plt.xlim(30,300)
    ##plt.xscale('log')
    #plt.ylim(40,80)

    ax1.set_title(title, fontsize=14)
    ax1.grid()
    ax1.legend(fontsize=14,loc=4)

    #plt.title(title, fontsize=14)
    #plt.grid()
    #plt.legend(fontsize=14)

    if outfile:
        plt.savefig(outfile)

    plt.close()

    return

if __name__ == '__main__':
    parser = OptionParser(__doc__)
    parser.add_option('--distance', default=None, help="LVC luminosity distance in Mpc")
    parser.add_option('--distance_err', default=None, help="LVC luminosity distance standard deviation in Mpc")
    parser.add_option('--time_delay',  default=None, help="Time since the merger in hours")
    parser.add_option('--fraction', default=90, help="Fraction of models you want to detect")
    parser.add_option('--magplot_file', default=None, help="Outfile for mag plot")
    parser.add_option('--expplot_file', default=None, help="Outfile for exp plot")
    parser.add_option('--report_file', default=None, help="Prefix for reports")
    parser.add_option('--filter', default=None, help="Single band for mag calculation")
    options, args = parser.parse_args(sys.argv[1:])

    if not options.distance:
        print("ERROR: You must specify the distance in Mpc with the --distance flag")
        sys.exit()

    if not options.distance_err:
        print("ERROR: You must specify the distance std in Mpc with the --distance_err flag")
        sys.exit()

    if not options.time_delay:
        print("ERROR: You must specify the time since merger in hours with the --time_delay flag")
        sys.exit()

    if options.report_file:
        suffix = os.path.splitext(options.report_file)[1]
        if (suffix == ".txt" or suffix == ".csv" ) :
            print("WARNING: report file argument is not supposed to include a .txt, .csv")
            
    if options.filter:
        if options.filter not in ['g', 'r', 'i', 'z']:
            print("ERROR: filter argument must be in ['g', 'r', 'i', 'z']")
            sys.exit()

    kn_calc = KNCalc(float(options.distance), float(options.distance_err), float(options.time_delay))

    if options.report_file:
        blue, red = gw170817(kn_calc.template_df_full)
        print_dict(blue, "GW170817-blue", outfile=options.report_file)
        print_dict(red, "GW170817-red", outfile=options.report_file)


    percentile_dict = calc_mag_fractions(kn_calc.template_df_full)
    cutoffs = mags_of_percentile(float(options.fraction), percentile_dict)
    cutoff_dict = {'%s_mag' %k : v for k, v in cutoffs.items()}
    for band in ['g', 'r', 'i', 'z']:
        cutoff_dict['%s_magerr' %band] = 0.00
    
    if options.filter:
        make_output_csv(np.linspace(0., 100., 101), percentile_dict, outfile=options.report_file, write_answer=True, flt=options.filter, fraction=options.fraction)
    else:
        make_output_csv(np.linspace(0., 100., 101), percentile_dict, outfile=options.report_file)
    
    if float(options.fraction) < 1.0:
        print_dict(cutoff_dict, "%.2f Detection Probability Magnitude Thresholds" %float(options.fraction * 100), outfile=options.report_file)
    else:
        print_dict(cutoff_dict, "%.2f Detection Probability Magnitude Thresholds" %float(options.fraction), outfile=options.report_file)

    plot_title = "%s +/- %s Mpc  -- %.2f Days After Merger" %(options.distance, options.distance_err, float(options.time_delay) / 24.0)
    if options.magplot_file:
        make_plot(percentile_dict, blue, red, title=plot_title, outfile=options.magplot_file, fraction=options.fraction)
    if options.expplot_file:
        make_exptime_plot(percentile_dict, title=plot_title, outfile=options.expplot_file)





