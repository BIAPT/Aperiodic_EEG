#!/usr/bin/env python

from antropy import lziv_complexity
import pandas as pd
import numpy as np
import argparse

import mne
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate the aperiodic signal portion using different Models."
    )
    parser.add_argument(
        "data_dir", type=str, action="store", help="folder name containing the data in .fif format"
    )
    parser.add_argument(
        "output_dir", type=str, action="store", help="directory for results to be saved"
    )
    parser.add_argument(
        "patient_information",
        type=str,
        action="store",
        help="path to txt with information about participants",
    )
    parser.add_argument(
        "condition",
        type=str,
        action="store",
        help='The "task" or conditions you want to caluclate for example Baseline or Anesthesia',
    )

    # read out arguments
    args = parser.parse_args()
    cond = args.condition
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    output = []
    output_space = []

    # load patient info
    info = pd.read_csv(args.patient_information, sep="\t")
    P_IDS = info["Patient"]

    for P_ID in P_IDS:
        # print task number
        print("######################################## PARTICIPANT : ", P_ID)

        """
        1) LOAD Participant DATA
        """
        # load the epoched raw data
        input_fname = f"{args.data_dir}/sub-{P_ID}/eeg/epochs_{P_ID}_{cond}.fif"
        raw = mne.read_epochs(input_fname)
        # remove channels marked as bad and non-brain channels
        raw.drop_channels(raw.info["bads"])

        raw_filtered = raw.filter(1, 45)

        # define empty DataFrames to save Baseline data
        LZC = []

        if len(raw_filtered) > 30:
            epochs_cropped = raw_filtered[-30:]
        else:
            epochs_cropped = raw_filtered.copy()

        # loop over every time segment
        for i, w in enumerate(epochs_cropped):
            median = np.median(w, axis=1)
            # initiate dataframe to collect data for one window
            LZC_window = []
            for e in range(len(median)):

                bin_sig = (w[e] > median[e]).astype(np.int_)
                LZC_window.append(lziv_complexity(bin_sig, normalize=True))
            LZC.append(LZC_window)
            print(f"LZC calculated for window {i + 1} / {len(epochs_cropped)}")

        print(f"Finished Subject  {P_ID}")

        LZC = np.array(LZC)

        output.append([P_ID, np.mean(LZC, axis=(0, 1))])
        output_space.append([P_ID, np.mean(LZC, axis=0)])

    output_df = pd.DataFrame(output, columns=["ID", f"LZC_{cond}"])
    output_space_df = pd.DataFrame(output_space, columns=["ID", f"LZC_{cond}"])

    output_space_df.to_csv(f"{output_dir}/Complexity_space_{cond}.txt", index=False, sep=" ")
    output_df.to_csv(f"{output_dir}/Complexity_{cond}.txt", index=False, sep=" ")
