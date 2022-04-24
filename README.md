# One_over_f
 Project by Charlotte Maschke



PART 1:

The first part contains the Analysis of data from the Baseline recording.

1. Calculate the PSD
'python _1_spectral_decomposition.py data_aperiodic RESULTS/spectrum_Baseonly data_aperiodic/data_baseonly.txt Base'
This outputs:
- Frequency.txt   > containing the frequencies for every PSD, depending on the Mehod
- PSD.pkl         > contains all PSDS for every epoch, patient and channel
- summary.pdf     > Summary of patients individual PSD

2. Parametrize PSD
'python _2a_parametrization.py RESULTS/spectrum_Baseonly RESULTS/aperiodic_Baseonly data_aperiodic data_aperiodic/data_baseonly.txt Base'
This outputs:
- Parametrized_Base.txt   > contains all parameters of interest for this paper
- PSDS_aper_Base.txt      > The Aperiodic part of the PSD for every subject
- PSDS_flat_Base.txt      > The Oscillation part of the PSD for every subject
- PSDS_norm_Base.txt      > The 'Normal' PSD for every subject
- space                   > contains all model parameters in a space-resolved way (will be used by topological plot)
- summary_param.pdf       > Individual patients model Fit Plots

2. (optional) Step 2 can also be done time-resolved (This takes MUCH more time)
'python _2a_parametrization.py RESULTS/spectrum_Baseonly RESULTS/aperiodic_Baseonly data_aperiodic data_aperiodic/data_baseonly.txt Base'
This outputs:
- Parametrized_Base.txt   > contains all parameters of interest for this paper
- PSDS_aper_Base.txt      > The Aperiodic part of the PSD for every subject
- PSDS_flat_Base.txt      > The Oscillation part of the PSD for every subject
- PSDS_norm_Base.txt      > The 'Normal' PSD for every subject
- space                   > contains all model parameters in a space-resolved way (will be used by topological plot)
- summary_param.pdf       > Individual patients model Fit Plots


3. Calculate Bandpower:
'python _2b_bandpower.py RESULTS/aperiodic_Baseonly data_aperiodic/data_baseonly.txt Base'
This outputs:
- Band_Power_Base.txt     > which contains the absolute Bandpower per participant

4. Calculate Complexity:
'_2c_complexity.py data_aperiodic RESULTS/aperiodic_baseonly/Complexity data_aperiodic/data_baseonly.txt Base'
This outputs:
- Complexity.txt          > the complexity for every participant
- Complexity_space.txt    > the complexity for every participant space-resolved for later plotting



PART2: Run the same thing for 2 staes (Base and Anes)
python _1_spectral_decomposition.py data_aperiodic RESULTS/spectrum_2states data_aperiodic/data_2states.txt Base
python _1_spectral_decomposition.py data_aperiodic RESULTS/spectrum_2states data_aperiodic/data_2states.txt Anes

python _2a_parametrization.py RESULTS/spectrum_2states RESULTS/aperiodic_2states data_aperiodic data_aperiodic/data_2states.txt Base
python _2a_parametrization.py RESULTS/spectrum_2states RESULTS/aperiodic_2states data_aperiodic data_aperiodic/data_2states.txt Anes

python _2b_bandpower.py RESULTS/aperiodic_2states data_aperiodic/data_2states.txt Base   
python _2b_bandpower.py RESULTS/aperiodic_2states data_aperiodic/data_2states.txt Anes

python _2c_complexity.py data_aperiodic RESULTS/aperiodic_2states/Complexity data_aperiodic/data_2states.txt Base
python _2c_complexity.py data_aperiodic RESULTS/aperiodic_2states/Complexity data_aperiodic/data_2states.txt Anes

python _3_topographic_plot.py data_aperiodic RESULTS/aperiodic_2states data_aperiodic/data_2states.txt -cond Base Anes -freq 30 45
python _3_topographic_plot.py data_aperiodic RESULTS/aperiodic_2states data_aperiodic/data_2states.txt -cond Anes Anes -freq 30 45
