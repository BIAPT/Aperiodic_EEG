#!/usr/bin/env python

from mne.time_frequency import psd_multitaper, psd_welch
import matplotlib.backends.backend_pdf as pltpdf
from utils.visualize import plot_two_curves
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import pickle
import mne
import os


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculates the Power spectral density.")
    parser.add_argument(
        "input_dir",
        type=str,
        action="store",
        help="folder name containing the data in epoched .fif data in BIDS format",
    )
    parser.add_argument(
        "output_dir", type=str, action="store", help="folder name where to save the power spectra"
    )
    parser.add_argument(
        "patient_information",
        type=str,
        action="store",
        help="path to txt with information about participants",
    )
    parser.add_argument(
        "condition", type=str, action="store", help='The "task" or condition you want to analyze'
    )
    args = parser.parse_args()

    """
           1)    PREPARE IN-AND OUTPUT
    """

    # make ouput directory
    output_dir = os.path.join(args.output_dir, args.condition)
    os.makedirs(output_dir, exist_ok=True)

    # prepare output pdf
    pdf = pltpdf.PdfPages(f"{output_dir}/spectral_decomposition_{args.condition}.pdf")

    # load patient IDS
    info = pd.read_csv(args.patient_information, sep="\t")
    P_IDS = info["Patient"]

    # define empty lists to save
    power_spectra_welch = []  # 1) Welch Power density
    power_spectra_mt = []  #    2) Multitaper power density

    for p_id in P_IDS:
        """
        2)      IMPORT DATA
        """
        # define epoch name
        input_fname = f"{args.input_dir}/sub-{p_id}/eeg/epochs_{p_id}_{args.condition}.fif"
        raw_epochs = mne.read_epochs(input_fname)

        # remove channels marked as bad and non-brain channels
        raw_epochs.drop_channels(raw_epochs.info["bads"])

        # crop data if necessary
        if len(raw_epochs) > 30:
            epochs_cropped = raw_epochs[-30:]
        else:
            epochs_cropped = raw_epochs.copy()

        """
            3)    Plot PSD for all electrodes
        """
        # show a general PSD
        fig, ax = plt.subplots(nrows=1, ncols=1)
        epochs_cropped.plot_psd(ax=ax, fmin=0.5, fmax=50, dB=True, show=False, average=False)
        ax.set_title(p_id + "  PSD")
        pdf.savefig(fig)
        plt.close()

        """
            4)    Compute PSD with Welch Method
        """
        # 1) using the Welch Method
        # Hanning window of 2 s
        # 1 s overlap
        sampling_rate = 250
        psds_welch, freqs_welch = psd_welch(
            epochs_cropped,
            fmin=0.5,
            fmax=50,
            n_fft=2 * sampling_rate,
            n_per_seg=2 * sampling_rate,
            n_overlap=sampling_rate,
        )
        psds_welch_db = 10 * np.log10(psds_welch)  # convert to dB

        # save the PSD for later analysis
        power_spectra_welch.append(psds_welch)

        """
            5)    Compute PSD with Multitaper Method
        """

        # 1) using the Multitaper Method
        # bandwith (frequency smoothing) of 0.5
        # resulting number of tapers: 4
        # See this paper: https://journals.physiology.org/doi/full/10.1152/physiol.00062.2015

        psds_mt, freqs_mt = psd_multitaper(epochs_cropped, fmin=0.5, fmax=50, bandwidth=1)
        psds_mt_db = 10 * np.log10(psds_mt)  # convert to dB

        # save the PSD for later analysis
        power_spectra_mt.append(psds_mt)

        """
            6)    Plot PSD Welch and Multitaper
        """
        # plot results
        fig = plot_two_curves(
            x1=freqs_welch,
            x2=freqs_mt,
            y1=psds_welch_db.mean(0),
            y2=psds_mt_db.mean(0),
            c1="green",
            c2="orange",
            l1="Welch",
            l2="Multitaper",
            title=f"{p_id} Power spectral density",
            lx="Frequency (Hz)",
            ly="Power Spectral Density (dB)",
        )
        pdf.savefig(fig)
        plt.close(fig)
        plt.show()

    # Add IDs at the end of the list
    power_spectra_mt.append(P_IDS)
    power_spectra_welch.append(P_IDS)

    with open(f"{output_dir}/PSD_Welch_{args.condition}.pkl", "wb") as f:
        pickle.dump(power_spectra_welch, f)

    with open(f"{output_dir}/PSD_Multitaper_{args.condition}.pkl", "wb") as f:
        pickle.dump(power_spectra_mt, f)

    freqs_mt = pd.DataFrame(freqs_mt)
    freqs_welch = pd.DataFrame(freqs_welch)
    freqs_welch.to_csv(
        "{}/Frequency_Welch_{}.txt".format(output_dir, args.condition),
        index=None,
        header=None,
        sep=" ",
    )
    freqs_mt.to_csv(
        "{}/Frequency_Multitaper_{}.txt".format(output_dir, args.condition),
        index=None,
        header=None,
        sep=" ",
    )

    pdf.close()
