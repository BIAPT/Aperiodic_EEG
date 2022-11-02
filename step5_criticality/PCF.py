#!/usr/bin/env python
from scipy.signal import hilbert
import pandas as pd
import numpy as np
import argparse
import mne.io
import mne

import os

def pcf(data):
    """Estimate the pair correlation function (PCF) in a network of
    oscillators, equivalent to the susceptibility in statistical physics.
    The PCF shows a sharp peak at a critical value of the coupling between
    oscillators, signaling the emergence of long-range correlations between
    oscillators.

    Parameters
    ----------
    data : 2d array
        The filtered input data, where the first dimension is the different
        oscillators, and the second dimension is time.

    Returns
    -------
    pcf : float
        The pair correlation function, a scalar value >= 0.
    orpa: float
        Absolute value of the order parameter (degree of synchronicity)
        averaged over time, being equal to 0 when the oscillators’ phases are
        uniformly distributed in [0,2π ) and 1 when they all have the same
        phase.
    orph_vector: 1d array (length = N_timepoints)
        Order parameter phase for every moment in time.
    orpa_vector: 1d array (length = N_timepoints)
        Absolute value of the order parameter (degree of synchronicity) for
        every moment in time.

    References
    ----------
    Yoon et al. (2015) Phys Rev E 91(3), 032814.
    """
    N_ch = data.shape[0]  # Nr of channels

    # calculate Phase of signal
    inst_phase = np.angle(hilbert(data))
    # get global synchronization order parameter z over time
    z_vector = np.mean(np.exp(1j * inst_phase), axis = 0)
    # get order phases
    orph_vector = np.arctan2(np.imag(z_vector), np.real(z_vector))
    #  r =|z| degree of synchronicity
    orpa_vector = np.abs(z_vector)

    # get PCF = variance of real part of order parameter
    # var(real(x)) == (mean(square(real(x))) - square(mean(real(x))))
    pcf = N_ch * np.var(np.real(z_vector))
    # time-averaged Order Parameter
    orpa = np.mean(orpa_vector);

    return pcf, orpa, orph_vector, orpa_vector


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate Criticality using different methods')
    parser.add_argument('data_dir', type=str, action='store',
                        help='folder name containing the data in .fif format')
    parser.add_argument('output_dir', type=str, action='store',
                        help='directory for results to be saved')
    parser.add_argument('part_info', type=str, action='store',
                        help='path to txt with information about participants')
    parser.add_argument('conditions', type=str, action='store',
                        help='path to txt with all condition names')


    # read out arguments
    args = parser.parse_args()
    output_dir = args.output_dir

    # make output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # load patient info and conditions
    info = pd.read_csv(args.part_info,sep = '\t', index_col=None)
    conds = pd.read_csv(args.conditions,sep = '\t') ['conditions']
    P_IDS = info['Patient']

    #loop over all conditions and particiants
    for c_tmp in conds:

        PCF_mean = []
        PCF_std = []
        OR_mean = []
        OR_std = []
        ID = []


        for p_id in P_IDS:
            print("Analyzing Criticality of {} in {}", p_id, c_tmp);

            ID_PCF = []
            ID_OR = []

            #################################
            #          LOAD  DATA          #
            #################################

            # load the epoched raw data
            input_fname = "{}/sub-{}/eeg/epochs_{}_{}.fif".format(args.data_dir, p_id, p_id, c_tmp)
            epochs = mne.read_epochs(input_fname)

            epochs.load_data()

            # prepare data
            epochs = epochs.filter(7, 13)

            # split data into non-overlaping 10 s windows
            nr_trials = min([len(epochs),30]);
            nr_channels =  epochs.info['nchan']

            # loop over every time segment
            for trial in range(nr_trials):
                #################################
                #          Calculate PCF        #
                #################################

                data_trial = epochs[trial].get_data()[0]

                # data needs to be channel x time
                if data_trial.shape[0] > data_trial.shape[1]:
                    data_trial=data_trial.T

                pcf_tmp, orpa_tmp, _, _ = pcf(data_trial)

                ID_PCF.append(pcf_tmp)
                ID_OR.append(orpa_tmp)

                print('Done {} Trial {}'.format(p_id, str(trial)))

            # fill in values
            PCF_mean.append(np.mean(ID_PCF))
            PCF_std.append(np.std(ID_PCF))
            OR_mean.append(np.mean(ID_OR))
            OR_std.append(np.std(ID_OR))

            ID.append(p_id);

        # save dataframe
        output_df = {'ID':ID, 'PCF_mean':PCF_mean, 'OR_mean':OR_mean}
        output_df = pd.DataFrame(output_df)
        output_df.to_csv('{}/PCF_{}.txt'.format(output_dir, c_tmp), index=False, sep=',')
