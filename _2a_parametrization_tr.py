#!/usr/bin/env python

from fooof import FOOOFGroup, FOOOF, fit_fooof_3d
import pandas as pd
import numpy as np
import argparse
import pickle
import os


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Decomposes the signal using FOOOF in a time-resolved way"
    )
    parser.add_argument(
        "input_dir",
        type=str,
        action="store",
        help="folder name containing the PSD pickles from step 1",
    )
    parser.add_argument(
        "output_dir", type=str, action="store", help="directory for results to be saved"
    )
    parser.add_argument(
        "data_dir",
        type=str,
        action="store",
        help="folder name containing the data in .fif format in BIDS format",
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
    parser.add_argument(
        "--method",
        action="store",
        default="Multitaper",
        choices=("Multitaper", "Welch"),
        help="The method used for Spectral decomposition in Step 1",
    )

    # read out arguments
    args = parser.parse_args()
    method = args.method
    cond = args.condition

    output_dir = args.output_dir
    os.makedirs(args.output_dir, exist_ok=True)

    # load participant info
    info = pd.read_csv(args.patient_information, sep="\t")
    P_IDS = info["Patient"]

    # create empty output
    slope_1_45_knee = []
    offset_1_45_knee = []
    error_knee = []

    # load power spectral data
    PSD = pickle.load(open(f"{args.input_dir}/{cond}/PSD_{method}_{cond}.pkl", "rb"))
    PSD_ID = PSD[-1]

    # load frequencies
    datapath = os.path.join(args.input_dir, cond, f"Frequency_{method}_{cond}.txt")
    freqs = pd.read_csv(datapath, sep=" ", header=None)
    freqs = np.array(np.squeeze(freqs))

    for p_id in P_IDS:
        # select individual PSD values, depending on ID
        index_id = np.where(PSD_ID == p_id)[0][0]
        PSD_p = PSD[index_id]

        """
        Run FOOOF model in 3D
        """
        # power spectra across data epochs within subjects, as [n_epochs, n_channels, n_freqs]
        PSD_p_3d = np.array(PSD_p)

        offset_tmp = []
        exponent_tmp = []
        error_tmp = []

        # initiate FOOOF model:
        fg = FOOOFGroup(aperiodic_mode="knee", min_peak_height=0.1, max_n_peaks=10)
        # Fit the 3D array of power spectra
        fgs = fit_fooof_3d(fg, freqs, PSD_p_3d, freq_range=[1, 45], n_jobs=10)

        # collect results
        for fg in fgs:
            # Aperiodic parameters averaged over space
            exponent_tmp.append(np.nanmean(fg.get_params("aperiodic_params", "exponent")))
            offset_tmp.append(np.nanmean(fg.get_params("aperiodic_params", "offset")))
            error_tmp.append(np.nanmean(fg.get_params("error")))

        slope_1_45_knee.append(np.mean(exponent_tmp))
        offset_1_45_knee.append(np.mean(offset_tmp))
        error_knee.append(np.mean(error_tmp))

        # save data for one participant in a txt
        params_p = [[p_id, np.mean(exponent_tmp), np.mean(offset_tmp), np.mean(error_tmp)]]

        params_time_df = pd.DataFrame(
            params_p, columns=["ID", "slope_1_45_knee", "offset_1_45_knee", "error_1_45_knee"]
        )

        params_time_df.to_csv(f"{output_dir}/Params_time_{p_id}_{cond}.txt", index=False, sep=" ")
        print(f"Finished Subject  {p_id}")

    params_df = pd.DataFrame()
    params_df["ID"] = P_IDS

    params_df["slope_1_45_knee"] = slope_1_45_knee
    params_df["offset_1_45_knee"] = offset_1_45_knee
    params_df["error_knee"] = error_knee

    params_df.to_csv(f"{output_dir}/Parametrized_{cond}_TR.txt", index=False, sep=" ")
