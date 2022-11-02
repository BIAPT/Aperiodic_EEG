### This repository contains the code for the paper: <strong> Aperiodic brain activity and response to anesthesia vary in disorders of consciousness </strong>

Project by Charlotte Maschke

Please check out the [preprint](https://www.biorxiv.org/content/10.1101/2022.04.22.489199v2) for this work

### Short description
Canonical analysis of human brain activity through electroencephalography (EEG) separates periodic and aperiodic components of the signal: the periodic component is further analyzed, while the aperiodic component is traditionally discarded. In this work, we show that the aperiodic component of EEG contains important information about levels of human consciousness.  Using EEG recorded from humans in pathological states of unconsciousness, we show that the aperiodic EEG component’s response to anesthesia varies with the individual’s level of consciousness. We further demonstrate that alterations in the aperiodic EEG component reflects the brain’s loss of network criticality and complexity, providing further evidence to support emerging theories that link criticality to mechanisms underpinning human consciousness.

Below, you find all code used to produce the results of this paper.

The criticality and complexity analysis code was originally provided by [Toker et al. (2022)](https://www.pnas.org/doi/abs/10.1073/pnas.2024455119). Only minor changes have been made to the complexity code in order to use it with two normalization methods.


The first part contains the Analysis of data from the Baseline recording only.

#### Setup

Install the required packages using `pip install -r requirements.txt`

#### PART 1: Calculate the PSD

`python _1_spectral_decomposition.py`
Inputs: please type `_1_spectral_decomposition.py --help` to get details on the required input

This outputs:
- Frequency.txt   > containing the frequencies for every PSD, depending on the method
- PSD.pkl         > contains all PSDS for every epoch, patient and channel
- summary.pdf     > Summary of patients individual PSD

#### PART 2: Parametrize PSD

`python _2a_parametrization.py`
Inputs: please type `_2a_parametrization.py --help` to get details on the required input

This outputs:
- Parametrized_Base.txt   > contains all parameters of interest for this paper
- PSDS_aper_Base.txt      > The Aperiodic part of the PSD for every subject
- PSDS_flat_Base.txt      > The Oscillation part of the PSD for every subject
- PSDS_norm_Base.txt      > The 'Normal' PSD for every subject
- space                   > contains all model parameters in a space-resolved way (will be used by topological plot)
- summary_param.pdf       > Individual patients model Fit Plots

(optional) Step 2 can also be done time-resolved (This takes MUCH more time).
Inputs: please type `python _2a_parametrization_tr.py --help` to get details on the required input
The output is the same as above

#### PART 3: Calculate Bandpower:
Inputs: please type `python _2b_bandpower.py --help` to get details on the required input

This outputs:
- Band_Power_Base.txt     > which contains the absolute Bandpower per participant

#### PART 4: Calculate Complexity:
The location step4_complexity contains a Matlab script `Calculate_all_complexity.m` it calls the complexity code which was provided by [Toker et al. (2022)](https://www.pnas.org/doi/abs/10.1073/pnas.2024455119)

#### PART 5: Calculate Criticality:
The location step5_criticality contains a Matlab script `Calculate_01criticality.m` it calls the modified 01 chaos test code which was provided by [Toker et al. (2022)](https://www.pnas.org/doi/abs/10.1073/pnas.2024455119)
The PCF can be calculated using `PCF.py`. Please type `python PCF.py --help` to get details on the required input. Both scripts output text files with the corresponding measure.
