%% Charlotte Maschke June 2022
% This script goal is to calculate complexity measures with the code published by Toker et.al 
% https://figshare.com/articles/software/Consciousness_is_supported_by_near-critical_cortical_electrodynamics/12949355

% Source Setup
%
INPUT_DIR = 'Users/charlotte/Documents/GitHub/one_over_f/data_aperiodic/';
OUTPUT_DIR = 'RESULTS/Complexity/';
CONDITION = 'Anes';
PART_FILE = "data_aperiodic/data_2states.txt";

% load participant info
opts = detectImportOptions(PART_FILE,'Delimiter','\t');
info = readtable(PART_FILE,opts);
P_IDS = info.Patient;

univ_phasen_LZC = {};
univ_shufn_LZC = {};
conc_phasen_LZC = {};
conc_shufn_LZC = {};
ID = {};

%% loop over all particiopants and stepsizes and calculate LZC
for p = 1:length(P_IDS)
    p_id = P_IDS{p};
    
    fprintf("Analyzing complexity of '%s' in '%s' \n", p_id,CONDITION);

    participant_in= strcat('sub-', p_id, '_task-',CONDITION,'_eeg.set');
    part_dir = strcat(INPUT_DIR,'sub-', p_id, filesep, 'eeg', filesep);

    %% Load data
    EEG = pop_loadset('filename',participant_in,'filepath',part_dir);

    %recording = load_set(participant_in, 'data/raw/');
    EEG = pop_resample( EEG, 250);
  
    fs = EEG.srate;

    % Filter the data Low pass 45 Hz
    lowpass = 45;
    EEG = pop_eegfiltnew(EEG, [], lowpass, [], false, [], 0); % Lowpass filter
  
    data = EEG.data;

    % split data into non-overlaping 10 s windows
    window_size = 10; % in seconds
    step_size = 10; % in seconds
    [Epochs] = create_window(data, fs, window_size, step_size);
    
    % make temporal Collection over all trials
    tmp_univ_shuf = {};
    tmp_univ_phase = {};
    tmp_conc_shuf = {};
    tmp_conc_phase = {};

    % take same nr of trials, maximal 30
    nr_trials = min([size(Epochs,1),30]);

    parfor trial=1:nr_trials
        % Schartner method for Conc LZC with 2 types of Normalization
        data_trial = squeeze(Epochs(trial,:,:));
        [pn_LZC_concat, ] = fJLZC(data_trial,fs,'concat','phase-rand');
        [sn_LZC_concat, ] = fJLZC(data_trial,fs,'concat','shuffle');
        tmp_conc_phase = [tmp_conc_phase, pn_LZC_concat];
        tmp_conc_shuf = [tmp_conc_shuf, sn_LZC_concat];

        % Univariate Method
        data_trial = squeeze(Epochs(trial,:,:));
        tmp_ch_pn = {}; %phase norm
        tmp_ch_sn = {}; %shuffle norm
        for ch =1:size(data,1)
            [pn_LZC_univ, ] = fJLZC(data_trial(ch,:),fs, '','phase-rand');
            [sn_LZC_univ, ] = fJLZC(data_trial(ch,:),fs, '','shuffle');
            tmp_ch_pn = [tmp_ch_pn, pn_LZC_univ];
            tmp_ch_sn = [tmp_ch_sn, sn_LZC_univ];
        end

        tmp_univ_phase = [tmp_univ_phase, median(cell2mat(tmp_ch_pn))];
        tmp_univ_shuf = [tmp_univ_shuf, median(cell2mat(tmp_ch_sn))];
        
        display("Finished Baseline " + string(p_id) + "trial " +string(trial))
    end
    
    
    % Fill in values to save
    univ_phasen_LZC = [univ_phasen_LZC, mean(cell2mat(tmp_univ_phase))];
    univ_shufn_LZC  = [univ_shufn_LZC, mean(cell2mat(tmp_univ_shuf))];
    conc_phasen_LZC = [conc_phasen_LZC, mean(cell2mat(tmp_conc_phase))];
    conc_shufn_LZC  = [conc_shufn_LZC, mean(cell2mat(tmp_conc_shuf))];
    ID = [ID, p_id];

    %% save data
    T = table(ID(:), univ_phasen_LZC(:), univ_shufn_LZC(:), conc_phasen_LZC(:), conc_shufn_LZC(:),...
        'VariableNames', { 'ID', 'univ_phasen_LZC', 'univ_shufn_LZC', 'conc_phasen_LZC', 'conc_shufn_LZC'});
    % Write data to text file
    writetable(T, strcat(OUTPUT_DIR,'univ_complexity_',CONDITION,'.txt'))
end

%% save data
T = table(ID(:), univ_phasen_LZC(:), univ_shufn_LZC(:), conc_phasen_LZC(:), conc_shufn_LZC(:),...
    'VariableNames', { 'ID', 'univ_phasen_LZC', 'univ_shufn_LZC', 'conc_phasen_LZC', 'conc_shufn_LZC'});
% Write data to text file
writetable(T, strcat(OUTPUT_DIR,'univ_complexity_',CONDITION,'.txt'))
