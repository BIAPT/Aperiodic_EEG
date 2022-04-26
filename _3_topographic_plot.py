#!/usr/bin/env python

import matplotlib.backends.backend_pdf as pltpdf
from electrode_list import electrode_list
import matplotlib.pyplot as plt
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
        "data_dir",
        type=str,
        action="store",
        help="folder name containing the data in .fdt and .set format",
    )
    parser.add_argument(
        "input_dir", type=str, action="store", help="folder name containing the data to plot"
    )
    parser.add_argument(
        "patient_information",
        type=str,
        action="store",
        help="path to txt with information about participants",
    )
    parser.add_argument(
        "--conditions",
        "-cond",
        nargs="*",
        action="store",
        default="Baseline Anesthesia",
        help='The "task" or conditions you want to compare for example Baseline Anesthesia'
        "can be only Base or Baseline and Anesthesia",
    )
    parser.add_argument(
        "--frequency_range",
        "-freq",
        nargs="*",
        action="store",
        default="1 40",
        help="The freqency band to calculate the aperiodic signal on. For example 1 20",
    )
    parser.add_argument(
        "--method",
        action="store",
        default=["Multitaper"],
        choices=("Multitaper", "Welch"),
        help="The method used for Spectral decomposition in Step 1",
    )

    args = parser.parse_args()
    nr_cond = len(args.conditions)
    frequency_range = [int(args.frequency_range[0]), int(args.frequency_range[1])]
    input_dir = args.input_dir
    method = args.method

    # make ouput directory
    output_dir = os.path.join(input_dir, "topoplot")
    os.makedirs(output_dir, exist_ok=True)

    # prepare output pdf
    pdf = pltpdf.PdfPages(
        f"{output_dir}/topoplot_aperiodic_signal_{frequency_range[0]}_{frequency_range[1]}.pdf"
    )

    col_names = electrode_list

    # load patient info
    info = pd.read_csv(args.patient_information, sep="\t")
    P_IDS = info["Patient"]

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
        params_B = pd.read_csv(
            f"{input_dir}/space/Params_space_{p_id}_{cond_B}_{frequency_range[0]}_{frequency_range[1]}.txt",
            sep=" ",
        )

        if nr_cond == 2:
            params_A = pd.read_csv(
                f"{input_dir}/space/Params_space_{p_id}_{cond_A}_{frequency_range[0]}_{frequency_range[1]}.txt",
                sep=" ",
            )

        # select p_id from aperiodic dataframe
        params_B_id = params_B[params_B["ID"] == p_id]
        # select only exponent data
        params_B_id = params_B_id["slope_space"].values[0]
        # remove the brakets from the string (some artifact from the saving format)
        params_B_id = params_B_id.replace("[", "")
        params_B_id = params_B_id.replace("]", "")
        exponent_B_id = params_B_id.split(sep=",")
        exponent_B_id = np.array(exponent_B_id).astype(float)
        exponent_B_id = exponent_B_id * -1

        if nr_cond == 2:
            # select p_id from aperiodic dataframe
            params_A_id = params_A[params_A["ID"] == p_id]
            # select only exponent data
            params_A_id = params_A_id["slope_space"].values[0]
            # remove the brakets from the string (some artifact from the saving format)
            params_A_id = params_A_id.replace("[", "")
            params_A_id = params_A_id.replace("]", "")
            exponent_A_id = params_A_id.split(sep=",")
            exponent_A_id = np.array(exponent_A_id).astype(float)
            exponent_A_id = exponent_A_id * -1

        # imput raw data (needed later for plotting only)
        input_fname = f"{args.data_dir}/sub-{p_id}/eeg/epochs_{p_id}_{cond_B}.fif"
        # remove channels marked as bad and non-brain channels
        raw_B = mne.read_epochs(input_fname)
        raw_B.drop_channels(raw_B.info["bads"])

        # select all channels to visualize
        exponent_id_select_B = list(exponent_B_id)

        if nr_cond == 2:
            input_fname = f"{args.data_dir}/sub-{p_id}/eeg/epochs_{p_id}_{cond_A}.fif"
            raw_A = mne.read_epochs(input_fname)
            # remove channels marked as bad and non-brain channels
            raw_A.drop_channels(raw_A.info["bads"])

            # compare electrodes to get overlying ones
            ch_A = np.array(raw_A.info.ch_names)
            ch_B = np.array(raw_B.info.ch_names)
            keep_A = np.isin(ch_A, ch_B)
            keep_B = np.isin(ch_B, ch_A)

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
        fig.suptitle(p_id + "__" + cond_B)
        pdf.savefig(fig)
        plt.close()

        # expand list to plot average at end
        for pos, el in enumerate(raw_B.info.ch_names):
            if col_names.__contains__(el):
                plot_all_Base[el][p_id] = exponent_id_select_B[pos]

        if nr_cond == 2:
            fig = plt.figure()
            im, _ = mne.viz.plot_topomap(
                exponent_id_select_A, raw_A.info, vmin=-4, vmax=0, show=False
            )
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
            exponent_diff = list(np.array(exponent_id_select_B) - np.array(exponent_id_select_A))

            fig = plt.figure()
            im, _ = mne.viz.plot_topomap(
                exponent_diff, raw_A.info, vmin=0, vmax=3, show=False, cmap="jet"
            )
            fig.colorbar(im)
            fig.subplots_adjust(top=0.8)
            fig.suptitle(p_id + "__" + "Base - Anes ")
            pdf.savefig(fig)
            plt.close()

            # expand list to plot average at end
            for pos, el in enumerate(raw_A.info.ch_names):
                if col_names.__contains__(el):
                    plot_all_diff[el][p_id] = exponent_diff[pos]

    plot_all_Base.to_csv(
        f"{output_dir}/Base_all_exp_{frequency_range[0]}_{frequency_range[1]}.txt",
        index=True,
        sep=";",
    )
    plot_all_Anes.to_csv(
        f"{output_dir}/Anes_all_exp_{frequency_range[0]}_{frequency_range[1]}.txt",
        index=True,
        sep=";",
    )
    plot_all_diff.to_csv(
        f"{output_dir}/Diff_all_exp_{frequency_range[0]}_{frequency_range[1]}.txt",
        index=True,
        sep=";",
    )

    pdf.close()
