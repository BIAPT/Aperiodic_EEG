#!/usr/bin/env python
from _2a_parametrization import get_power_band
import pandas as pd
import numpy as np
import argparse
import os


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculate the oscillatory part of the signal")
    parser.add_argument(
        "input_dir",
        type=str,
        action="store",
        help="folder name containing the spectrograms from step 2a",
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

    # make output directory
    output_dir = os.path.join(args.input_dir, "bandpower")
    os.makedirs(output_dir, exist_ok=True)

    # load patient info
    info = pd.read_csv(args.patient_information, sep="\t")
    P_IDS = info["Patient"]

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
    PSD = pd.read_csv(f"{args.input_dir}/PSDS_norm_{cond}.txt", sep=" ")
    PSD_flat = pd.read_csv(f"{args.input_dir}/PSDS_flat_{cond}.txt", sep=" ")
    # ignore less than 0 values
    PSD_flat[PSD_flat < 0] = 0

    # load frequencies
    freqs = np.arange(1, 45.1, 0.1)

    for i, p_id in enumerate(P_IDS):
        # select individual PSD values, depending on ID
        PSD_p = PSD.iloc[i, :]
        PSD_p_flat = PSD_flat.iloc[i, :]

        # delta
        power_delta.append(get_power_band(1, 4, freqs, PSD_p))
        flat_power_delta.append(get_power_band(1, 4, freqs, PSD_p_flat))

        # theta
        power_theta.append(get_power_band(4, 8, freqs, PSD_p))
        flat_power_theta.append(get_power_band(4, 8, freqs, PSD_p_flat))

        # alpha
        power_alpha.append(get_power_band(8, 13, freqs, PSD_p))
        flat_power_alpha.append(get_power_band(8, 13, freqs, PSD_p_flat))

        # beta
        power_beta.append(get_power_band(13, 30, freqs, PSD_p))
        flat_power_beta.append(get_power_band(13, 30, freqs, PSD_p_flat))

        # gamma
        power_gamma.append(get_power_band(30, 45, freqs, PSD_p))
        flat_power_gamma.append(get_power_band(30, 45, freqs, PSD_p_flat))

        print(f"Finished Subject  {p_id}")

    params_df = pd.DataFrame()
    params_df["ID"] = P_IDS
    params_df["power_delta"] = power_delta
    params_df["power_theta"] = power_theta
    params_df["power_alpha"] = power_alpha
    params_df["power_beta"] = power_beta
    params_df["power_gamma"] = power_gamma

    params_df["flat_power_delta"] = flat_power_delta
    params_df["flat_power_theta"] = flat_power_theta
    params_df["flat_power_alpha"] = flat_power_alpha
    params_df["flat_power_beta"] = flat_power_beta
    params_df["flat_power_gamma"] = flat_power_gamma

    params_df.to_csv(f"{output_dir}/Band_Power_{cond}.txt", index=False, sep=" ")
