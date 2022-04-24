#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse
import os


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate the oscillatory part of the signal')
    parser.add_argument('input_dir', type=str, action='store',
                        help='folder name containing the spectrograms from step 2a')
    parser.add_argument('patient_information', type=str, action='store',
                        help='path to txt with information about participants')
    parser.add_argument('condition', nargs=1, type=str, action='store',
                        help='The "task" or conditions you want to caluclate for example Baseline or Anesthesia')

    # read out arguments
    args = parser.parse_args()
    cond = args.condition[0]

    # make output directory
    output_dir = os.path.join(args.input_dir, 'bandpower')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # load patient info
    info = pd.read_csv(args.patient_information, sep='\t')
    P_IDS = info['Patient']

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

    # load power spectral data
    PSD = pd.read_csv("{}/PSDS_norm_{}.txt".format(args.input_dir, cond), sep = ' ')
    PSD_flat = pd.read_csv("{}/PSDS_flat_{}.txt".format(args.input_dir, cond), sep = ' ')
    # ignore less than 0 values
    PSD_flat[PSD_flat < 0] = 0

    # load frequencies
    freqs = np.arange(1,45.1, 0.1)

    for i, p_id  in enumerate(P_IDS):
        # select individual PSD values, depending on ID
        PSD_p = PSD.iloc[i,:]
        PSD_p_flat = PSD_flat.iloc[i,:]

        # delta
        index_toselect = np.where((freqs >= 1) & (freqs <= 4))[0]
        freqs_select = freqs[index_toselect]
        PSD_delta = PSD_p[index_toselect]
        PSD_delta_flat = PSD_p_flat[index_toselect]
        # calculate the absolute mean power in this frequency range:
        power_delta.append(np.mean(PSD_delta))
        flat_power_delta.append(np.mean(PSD_delta_flat))

        # theta
        index_toselect = np.where((freqs >= 4) & (freqs <= 8))[0]
        freqs_select = freqs[index_toselect]
        PSD_theta = PSD_p[index_toselect]
        PSD_theta_flat = PSD_p_flat[index_toselect]
        # calculate the absolute mean power in this frequency range:
        power_theta.append(np.mean(PSD_theta))
        flat_power_theta.append(np.mean(PSD_theta_flat))

        # alpha
        index_toselect = np.where((freqs >= 8) & (freqs <= 13))[0]
        freqs_select = freqs[index_toselect]
        PSD_alpha = PSD_p[index_toselect]
        PSD_alpha_flat = PSD_p_flat[index_toselect]
        # calculate the absolute mean power in this frequency range:
        power_alpha.append(np.mean(PSD_alpha))
        flat_power_alpha.append(np.mean(PSD_alpha_flat))

        # beta
        index_toselect = np.where((freqs >= 13) & (freqs <= 30))[0]
        freqs_select = freqs[index_toselect]
        PSD_beta = PSD_p[index_toselect]
        PSD_beta_flat = PSD_p_flat[index_toselect]
        # calculate the absolute mean power in this frequency range:
        power_beta.append(np.mean(PSD_beta))
        flat_power_beta.append(np.mean(PSD_beta_flat))

        # gamma
        index_toselect = np.where((freqs >= 30) & (freqs <= 45))[0]
        freqs_select = freqs[index_toselect]
        PSD_gamma = PSD_p[index_toselect]
        PSD_gamma_flat = PSD_p_flat[index_toselect]
        # calculate the absolute mean power in this frequency range:
        power_gamma.append(np.mean(PSD_gamma))
        flat_power_gamma.append(np.mean(PSD_gamma_flat))

        print("Finished Subject  {}".format(p_id))

    params_df = pd.DataFrame()
    params_df['ID'] = P_IDS
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

    params_df.to_csv('{}/Band_Power_{}.txt'.format(output_dir, cond), index=False, sep=' ')
