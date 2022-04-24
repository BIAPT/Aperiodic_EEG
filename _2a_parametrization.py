#!/usr/bin/env python

from sklearn.linear_model import LinearRegression
import matplotlib.backends.backend_pdf as pltpdf
import matplotlib.pyplot as plt
from fooof import FOOOFGroup, FOOOF
import pandas as pd
import numpy as np
import argparse
import pickle
import os

def extract_power(PSD_p, freqs, rel):

    # bring the PSD in positive area if it crosses the y-axis
    # as we only calculate the relative Power of every bandwidth, this does not matter
    if rel:
        if min(PSD_p) < 0:
            PSD_p = PSD_p + np.abs(min(PSD_p))

    # delta
    index_toselect = np.where((freqs >= 1) & (freqs <= 4))[0]
    PSD_delta = PSD_p[index_toselect]
    power_delta = np.mean(PSD_delta)

    # theta
    index_toselect = np.where((freqs >= 4) & (freqs <= 8))[0]
    PSD_theta = PSD_p[index_toselect]
    power_theta = np.mean(PSD_theta)

    # alpha
    index_toselect = np.where((freqs >= 8) & (freqs <= 13))[0]
    PSD_alpha = PSD_p[index_toselect]
    power_alpha = np.mean(PSD_alpha)

    # beta
    index_toselect = np.where((freqs >= 13) & (freqs <= 30))[0]
    PSD_beta = PSD_p[index_toselect]
    power_beta = np.mean(PSD_beta)

    # gamma
    index_toselect = np.where((freqs >= 30) & (freqs <= 45))[0]
    PSD_gamma = PSD_p[index_toselect]
    power_gamma = np.mean(PSD_gamma)

    power_total = power_delta + power_theta + power_alpha + power_delta + power_gamma

    if rel:
        power_delta = (power_delta/power_total)*100
        power_theta = (power_theta/power_total)*100
        power_alpha = (power_alpha/power_total)*100
        power_beta = (power_beta/power_total)*100
        power_gamma = (power_gamma/power_total)*100

    return power_delta, power_theta, power_alpha, power_beta, power_gamma


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate the oscillatory part of the signal')
    parser.add_argument('input_dir', type=str, action='store',
                        help='folder name containing the psd pickles from step 1')
    parser.add_argument('output_dir', type=str, action='store',
                        help='directory for results to be saved')
    parser.add_argument('data_dir', type=str, action='store',
                        help='folder name containing the data in .fif format')
    parser.add_argument('patient_information', type=str, action='store',
                        help='path to txt with information about participants')
    parser.add_argument('condition', nargs=1, type=str, action='store',
                        help='The "task" or conditions you want to caluclate for example Baseline or Anesthesia')
    parser.add_argument('--method', nargs=1, action='store', default=['Multitaper'], choices=('Multitaper','Welch'),
                        help='The method used for Spectral decomposition in Step 1')

    # read out arguments
    args = parser.parse_args()
    method = args.method[0]
    cond = args.condition[0]
    output_dir = args.output_dir

    # prepare output pdf
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    pdf = pltpdf.PdfPages("{}/summary_parametrization.pdf".format(output_dir))

    # load patient info
    info = pd.read_csv(args.patient_information, sep='\t')
    P_IDS = info['Patient']

    # define empty DataFrames values of interest
    peak_center = []
    peak_strength = []
    nr_peaks = []

    # define empty DataFrames to save Baseline data
    power_delta = []
    power_theta = []
    power_alpha = []
    power_beta = []
    power_gamma = []

    # define empty DataFrames to save Baseline data
    flat_power_delta = []
    flat_power_theta = []
    flat_power_alpha = []
    flat_power_beta = []
    flat_power_gamma = []

    slope_1_45_fixed = []
    offset_1_45_fixed = []
    r_squared_fixed = []
    error_fixed = []

    slope_30_45 = []
    offset_30_45 = []
    slope_30_45_linreg = []
    offset_30_45_linreg = []

    slope_1_45_knee = []
    offset_1_45_knee = []
    knee_1_45_knee = []
    r_squared_knee = []
    error_knee = []

    psds_norm = []
    psds_aper = []
    psds_flat = []


    # load power spectral data
    PSD = pickle.load(open("{}/{}/PSD_{}_{}.pkl".format(args.input_dir, cond, method, cond), "rb"))
    PSD_ID = PSD[-1]

    # load frequencies
    datapath = os.path.join(args.input_dir, cond, 'Frequency_{}_{}.txt'.format(method, cond))
    freqs = pd.read_csv(datapath, sep=' ', header=None)
    freqs = np.squeeze(freqs)

    for p_id in P_IDS:
        # select individual PSD values, depending on ID
        index_id = np.where(PSD_ID == p_id)[0][0]
        PSD_p = PSD[index_id]

        # average over time
        PSD_p = np.mean(PSD_p, axis=0)

        # load frequencies
        datapath = os.path.join(args.input_dir, cond, 'Frequency_{}_{}.txt'.format(method, cond))
        freqs = pd.read_csv(datapath, sep=' ', header=None)
        freqs = np.squeeze(freqs)
        freqs = np.array(freqs)

        """
        detekt peak in averaged PSD
        """
        # initialize fooof
        f = FOOOF(aperiodic_mode='knee', min_peak_height=0.5, max_n_peaks=1)
        f.fit(freqs, np.mean(PSD_p, axis=0), [1, 45])

        nr_peaks.append(f.n_peaks_)
        if f.n_peaks_ != 0:
            peak_center.append(f.get_params('peak_params', 'CF'))
            peak_strength.append(f.get_params('peak_params', 'PW'))
        else :
            peak_center.append(np.nan)
            peak_strength.append(np.nan)

        fig = f.plot()
        plt.title(p_id)
        pdf.savefig(fig)
        plt.close()

        """
        apply fooof object
        """
        slope_tmp = []
        offset_tmp = []
        knee_tmp = []
        error_tmp = []
        r2_tmp = []

        psd_norm = []
        psd_flat = []
        psd_aper = []

        # define empty DataFrames to save Baseline data
        power_delta_tmp = []
        power_theta_tmp = []
        power_alpha_tmp = []
        power_beta_tmp = []
        power_gamma_tmp = []

        # define empty DataFrames to save Baseline data
        flat_power_delta_tmp = []
        flat_power_theta_tmp = []
        flat_power_alpha_tmp = []
        flat_power_beta_tmp = []
        flat_power_gamma_tmp = []

        """
        Run FOOOF model in knee mode 1 -45 Hz for every electrode individually
        """
        for e in range(PSD_p.shape[0]):
            psd_tmp = PSD_p[e,:]

            # initialize fooof
            f = FOOOF(aperiodic_mode='knee', min_peak_height=0.1, max_n_peaks=10)
            f.fit(freqs, psd_tmp, [1, 45])

            # check whether model fit was sucessfull
            if not np.isnan(f.error_):
                slope_tmp.append(f.get_params('aperiodic_params', 'exponent'))
                offset_tmp.append(f.get_params('aperiodic_params', 'offset'))
                knee_tmp.append(f.get_params('aperiodic_params', 'knee'))
                error_tmp.append(f.get_params('error'))
                r2_tmp.append(f.get_params('r_squared'))

                psd_flat.append(f._spectrum_flat)
                psd_norm.append(f.power_spectrum)
                psd_aper.append(f._spectrum_peak_rm)

                pde, pt, pa, pb, pg = extract_power(f.power_spectrum, f.freqs, rel = True)
                power_delta_tmp.append(pde)
                power_theta_tmp.append(pt)
                power_alpha_tmp.append(pa)
                power_beta_tmp.append(pb)
                power_gamma_tmp.append(pg)

                pde, pt, pa, pb, pg = extract_power(f._spectrum_flat, f.freqs, rel = False)
                flat_power_delta_tmp.append(pde)
                flat_power_theta_tmp.append(pt)
                flat_power_alpha_tmp.append(pa)
                flat_power_beta_tmp.append(pb)
                flat_power_gamma_tmp.append(pg)

            else:
                slope_tmp.append(np.nan)
                offset_tmp.append(np.nan)
                knee_tmp.append(np.nan)

        # Fill in the values
        slope_1_45_knee.append(np.nanmean(slope_tmp))
        offset_1_45_knee.append(np.nanmean(offset_tmp))
        knee_1_45_knee.append(np.nanmean(knee_tmp))
        r_squared_knee.append(np.mean(r2_tmp))
        error_knee.append(np.mean(error_tmp))

        psds_norm.append(np.mean(psd_norm, axis =0))
        psds_aper.append(np.mean(psd_aper, axis =0))
        psds_flat.append(np.mean(psd_flat, axis =0))

        power_delta.append(np.mean(power_delta_tmp))
        power_theta.append(np.mean(power_theta_tmp))
        power_alpha.append(np.mean(power_alpha_tmp))
        power_beta.append(np.mean(power_beta_tmp))
        power_gamma.append(np.mean(power_gamma_tmp))

        flat_power_delta.append(np.mean(flat_power_delta_tmp))
        flat_power_theta.append(np.mean(flat_power_theta_tmp))
        flat_power_alpha.append(np.mean(flat_power_alpha_tmp))
        flat_power_beta.append(np.mean(flat_power_beta_tmp))
        flat_power_gamma.append(np.mean(flat_power_gamma_tmp))

        params_space = []
        params_space.append([p_id, offset_tmp, slope_tmp, knee_tmp])
        space_out_dir = '{}/space'.format(output_dir)

        params_space_df = pd.DataFrame(params_space, columns=['ID', 'offset_space', 'slope_space', 'knee_space'])

        if not os.path.exists(space_out_dir):
            os.makedirs(space_out_dir)

        params_space_df.to_csv('{}/space/Params_space_{}_{}_{}_{}.txt'.format(output_dir, p_id, cond, 1, 45),
                               index=False, sep=' ')

        """
        Fit FOOOF model in fixed mode (here we can use the fooofgroup)
        """
        # initialize fooof
        fg = FOOOFGroup(aperiodic_mode='fixed', min_peak_height=0.1, max_n_peaks=10)
        fg.fit(freqs, PSD_p, [1, 45])

        slope_1_45_fixed.append(np.nanmean(fg.get_params('aperiodic_params', 'exponent')))
        offset_1_45_fixed.append(np.nanmean(fg.get_params('aperiodic_params', 'offset')))
        error_fixed.append(np.nanmean(fg.get_params('error')))
        r_squared_fixed.append(np.nanmean(fg.get_params('r_squared')))

        """
        Fit FOOOF model for 30 to 45 Hz in fixed mode (here we can use the fooofgroup)
        """
        # initialize fooof
        fg = FOOOFGroup(aperiodic_mode='fixed', min_peak_height=0.1, max_n_peaks=10)
        fg.fit(freqs, PSD_p, [30, 45])

        slope_30_45.append(np.nanmean(fg.get_params('aperiodic_params', 'exponent')))
        offset_30_45.append(np.nanmean(fg.get_params('aperiodic_params', 'offset')))

        params_space = []
        params_space.append([p_id, list(fg.get_params('aperiodic_params', 'offset')),
                             list(fg.get_params('aperiodic_params', 'exponent'))])

        params_space_df = pd.DataFrame(params_space, columns=['ID', 'offset_space', 'slope_space'])

        params_space_df.to_csv('{}/space/Params_space_{}_{}_{}_{}.txt'.format(output_dir, p_id, cond, 30, 45),
                               index=False, sep=' ')

        """
            Fit the same 30-45 Hz Range with a linear regression instead
        """
        # create empty dataframe electrode x 1
        nr_chan = PSD_p.shape[0]
        exp_p = np.empty([nr_chan, 1])
        offs_p = np.empty([nr_chan, 1])

        # loop over electrodes to calculate linreg
        for e in range(nr_chan):
            # select only current channel
            PSD_tmp = PSD_p[e, :]

            # only use defined frequency band
            index_toselect = np.where((freqs >= 30) & (freqs <= 45))[0]
            freqs_select = freqs[index_toselect]
            PSD_tmp_select = PSD_tmp[index_toselect]

            # log the frequency and PSD
            freqs_select_log = np.array(np.log10(freqs_select))
            PSD_log = np.log10(PSD_tmp_select)

            #   Calculate Aperiodic signal Linear Regression Model
            mdl = LinearRegression()
            mdl.fit(freqs_select_log.reshape(-1, 1), PSD_log)

            # extract Aperiodic parameters
            offs_p[e] = mdl.intercept_
            exp_p[e] = mdl.coef_[0]

        offset_30_45_linreg.append(np.mean(offs_p))
        slope_30_45_linreg.append(np.mean(exp_p))

        print("Finished Subject  {}".format(p_id))



    params_df = pd.DataFrame()
    params_df['ID'] = P_IDS

    params_df['peak_freq'] = peak_center
    params_df['peak_power'] = peak_strength
    params_df['nr_peaks'] = nr_peaks

    params_df['slope_1_45_knee'] = slope_1_45_knee
    params_df['offset_1_45_knee'] = offset_1_45_knee
    params_df['knee_1_45_knee'] = knee_1_45_knee
    params_df['error_knee'] = error_knee
    params_df['r_squared_knee'] = r_squared_knee

    params_df['slope_1_45_fixed'] = slope_1_45_fixed
    params_df['offset_1_45_fixed'] = offset_1_45_fixed
    params_df['error_fixed'] = error_fixed
    params_df['r_squared_fixed'] = r_squared_fixed

    params_df['slope_30_45'] = slope_30_45
    params_df['offset_30_45'] = offset_30_45
    params_df['slope_30_45_linreg'] = slope_30_45_linreg
    params_df['offset_30_45_linreg'] = offset_30_45_linreg

    params_df['power_delta'] = power_delta
    params_df['power_theta'] = power_theta
    params_df['power_alpha'] = power_alpha
    params_df['power_beta'] = power_beta
    params_df['power_gamma'] = power_gamma

    params_df['flat_power_delta'] = flat_power_delta
    params_df['flat_power_theta'] = flat_power_theta
    params_df['flat_power_alpha'] = flat_power_alpha
    params_df['flat_power_beta'] = flat_power_beta
    params_df['flat_power_gamma'] = flat_power_gamma

    params_df.to_csv('{}/Parametrized_{}.txt'.format(output_dir, cond), index=False, sep=' ')

    pd.DataFrame(psds_norm).to_csv('{}/PSDS_norm_{}.txt'.format(output_dir, cond), index=False, sep=' ')
    pd.DataFrame(psds_aper).to_csv('{}/PSDS_aper_{}.txt'.format(output_dir, cond), index=False, sep=' ')
    pd.DataFrame(psds_flat).to_csv('{}/PSDS_flat_{}.txt'.format(output_dir, cond), index=False, sep=' ')

    pdf.close()
