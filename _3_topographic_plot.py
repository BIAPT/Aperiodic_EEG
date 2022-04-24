#!/usr/bin/env python

import matplotlib.backends.backend_pdf as pltpdf
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import mne
import os


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate the aperiodic signal portion using different Models.')
    parser.add_argument('data_dir', type=str, action='store',
                        help='folder name containing the data in .fdt and .set format')
    parser.add_argument('input_dir', type=str, action='store',
                        help='folder name containing the data to plot')
    parser.add_argument('patient_information', type=str, action='store',
                        help='path to txt with information about participants')
    parser.add_argument('--conditions', '-cond', nargs='*', action='store', default='Baseline Anesthesia',
                        help='The "task" or conditions you want to compare for example Baseline Anesthesia'
                             'can be only Base or Baseline and Anesthesia')
    parser.add_argument('--frequency_range', '-freq', nargs='*', action='store', default='1 40',
                        help='The freqency band to calculate the aperiodic signal on. For example 1 20')
    parser.add_argument('--method', nargs=1, action='store', default=['Multitaper'], choices=('Multitaper','Welch'),
                        help='The method used for Spectral decomposition in Step 1')

    args = parser.parse_args()
    nr_cond = len(args.conditions)
    frequency_range = [int(args.frequency_range[0]), int(args.frequency_range[1])]
    input_dir = args.input_dir
    method = args.method[0]

    # make ouput directory
    output_dir = os.path.join(input_dir, 'topoplot')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # prepare output pdf
    pdf = pltpdf.PdfPages("{}/topoplot_aperiodic_signal_{}_{}.pdf".format
                          (output_dir, frequency_range[0],frequency_range[1]))

    col_names = [ 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E9', 'E10', 'E11', 'E12', 'E13', 'E15',
                  'E16', 'E18', 'E19', 'E20', 'E22', 'E23', 'E24', 'E26', 'E27', 'E28', 'E29', 'E30',
                  'E31', 'E32', 'E33', 'E34', 'E35', 'E36', 'E37', 'E38', 'E39', 'E40', 'E41', 'E42',
                  'E44', 'E45', 'E46', 'E47', 'E50', 'E51', 'E52', 'E53', 'E54', 'E55', 'E57', 'E58',
                  'E59', 'E60', 'E61', 'E62', 'E64', 'E65', 'E66', 'E67', 'E69', 'E70', 'E71', 'E72',
                  'E74', 'E75', 'E76', 'E77', 'E78', 'E79', 'E80', 'E82', 'E83', 'E84', 'E85', 'E86',
                  'E87', 'E89', 'E90', 'E91', 'E92', 'E93', 'E95', 'E96', 'E97', 'E98', 'E100', 'E101',
                  'E102', 'E103', 'E104', 'E105', 'E106', 'E108', 'E109', 'E110', 'E111', 'E112', 'E114',
                  'E115', 'E116', 'E117', 'E118', 'E121', 'E122', 'E123', 'E124', 'E129']

    # load patient info
    info = pd.read_csv(args.patient_information,sep = '\t')
    P_IDS = info['Patient']

    # load empty frame to average at the end
    plot_all_Base = pd.DataFrame(columns=col_names, index=P_IDS)
    plot_all_Anes = pd.DataFrame(columns=col_names, index=P_IDS)
    plot_all_diff = pd.DataFrame(columns=col_names, index=P_IDS)

    # load conditions
    cond_B = str(args.conditions[0])
    if nr_cond == 2:
        cond_A = str(args.conditions[1])

    """
        1) load aperiodic data
    """
    for p_id in P_IDS:
        params_B = pd.read_csv('{}/space/Params_space_{}_{}_{}_{}.txt'.format
                               (input_dir, p_id, cond_B, frequency_range[0], frequency_range[1]),
                               sep=' ')

        if nr_cond == 2:
            params_A = pd.read_csv('{}/space/Params_space_{}_{}_{}_{}.txt'.format
                                   (input_dir, p_id, cond_A, frequency_range[0], frequency_range[1]),
                                   sep=' ')

        # select p_id from aperiodic dataframe
        params_B_id = params_B[params_B['ID']==p_id]
        # select only exponent data
        params_B_id = params_B_id['slope_space'].values[0]
        # remove the brakets from the string (some artifact from the saving format)
        params_B_id = params_B_id.replace('[', '')
        params_B_id = params_B_id.replace(']', '')
        exponent_B_id = params_B_id.split(sep=',')
        exponent_B_id = np.array(exponent_B_id).astype(float)
        exponent_B_id = exponent_B_id * -1


        if nr_cond == 2:
            # select p_id from aperiodic dataframe
            params_A_id = params_A[params_A['ID'] == p_id]
            # select only exponent data
            params_A_id = params_A_id['slope_space'].values[0]
            # remove the brakets from the string (some artifact from the saving format)
            params_A_id = params_A_id.replace('[', '')
            params_A_id = params_A_id.replace(']', '')
            exponent_A_id = params_A_id.split(sep=',')
            exponent_A_id = np.array(exponent_A_id).astype(float)
            exponent_A_id = exponent_A_id * -1

        # imput raw data (needed later for plotting only)
        input_fname = "{}/sub-{}/eeg/epochs_{}_{}.fif".format(args.data_dir, p_id, p_id, cond_B)
        # remove channels marked as bad and non-brain channels
        raw_B = mne.read_epochs(input_fname)
        raw_B.drop_channels(raw_B.info['bads'])

        # select all channels to visualize
        exponent_id_select_B = list(exponent_B_id)

        if nr_cond == 2:
            #input_fname = "{}/{}_{}.set".format(args.data_dir, p_id, cond_A)
            input_fname = "{}/sub-{}/eeg/epochs_{}_{}.fif".format(args.data_dir, p_id, p_id, cond_A)
            raw_A = mne.read_epochs(input_fname)
            # remove channels marked as bad and non-brain channels
            raw_A.drop_channels(raw_A.info['bads'])

            # compare electrodes to get overlying ones
            ch_A = np.array(raw_A.info.ch_names)
            ch_B = np.array(raw_B.info.ch_names)
            keep_A = np.isin(ch_A,ch_B)
            keep_B = np.isin(ch_B,ch_A)

            # select only the aperiodic components for both conditions
            exponent_id_select_B = list(exponent_B_id[keep_B])
            exponent_id_select_A = list(exponent_A_id[keep_A])

            # and do the same thing for the channel information
            raw_B.drop_channels(list(ch_B[np.invert(keep_B)]))
            raw_A.drop_channels(list(ch_A[np.invert(keep_A)]))

        fig = plt.figure()
        im, _ = mne.viz.plot_topomap(exponent_id_select_B, raw_B.info, vmin=-4, vmax=0, show=False)
        fig.colorbar(im)
        fig.subplots_adjust(top=0.8)
        fig.suptitle(p_id+"__"+cond_B)
        pdf.savefig(fig)
        plt.close()

        #expand list to plot average at end
        for pos, el in enumerate(raw_B.info.ch_names):
            if col_names.__contains__(el):
                plot_all_Base[el][p_id] = exponent_id_select_B[pos]

        if nr_cond == 2:
            fig = plt.figure()
            im, _ = mne.viz.plot_topomap(exponent_id_select_A, raw_A.info, vmin= -4, vmax=0, show=False)
            fig.colorbar(im)
            fig.subplots_adjust(top=0.8)
            fig.suptitle(p_id + "__" + cond_A)
            pdf.savefig(fig)
            plt.close()

            # expand list to plot average at end
            for pos, el in enumerate(raw_A.info.ch_names):
                if col_names.__contains__(el):
                    plot_all_Anes[el][p_id] = exponent_id_select_A[pos]

            # calculate the difference
            exponent_diff= list(np.array(exponent_id_select_B )-np.array(exponent_id_select_A))

            fig = plt.figure()
            im, _ = mne.viz.plot_topomap(exponent_diff, raw_A.info, vmin=0, vmax=3, show=False, cmap = 'jet')
            fig.colorbar(im)
            fig.subplots_adjust(top=0.8)
            fig.suptitle(p_id + "__" + 'Base - Anes ')
            pdf.savefig(fig)
            plt.close()

            # expand list to plot average at end
            for pos, el in enumerate(raw_A.info.ch_names):
                if col_names.__contains__(el):
                    plot_all_diff[el][p_id] = exponent_diff[pos]

    plot_all_Base.to_csv('{}/Base_all_exp_{}_{}.txt'.format(output_dir, frequency_range[0], frequency_range[1]),
                         index=True, sep=';')
    plot_all_Anes.to_csv('{}/Anes_all_exp_{}_{}.txt'.format(output_dir, frequency_range[0], frequency_range[1]),
                         index=True, sep=';')
    plot_all_diff.to_csv('{}/Diff_all_exp_{}_{}.txt'.format(output_dir, frequency_range[0], frequency_range[1]),
                         index=True, sep=';')

    pdf.close()

