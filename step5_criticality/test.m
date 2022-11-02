%% Charlotte Maschke June 2022
% This script goal is to calculate complexity measures with the code published by Toker et.al 
% https://figshare.com/articles/software/Consciousness_is_supported_by_near-critical_cortical_electrodynamics/12949355

% Source Setup
%
INPUT_DIR = 'Users/charlotte/Documents/GitHub/one_over_f/data_aperiodic/';
OUTPUT_DIR = 'RESULTS/complexity/';

% load participant info
opts = detectImportOptions("data_aperiodic/data_2states.txt",'Delimiter','\t');
info = readtable("data/raw/data_2states.txt",opts);
P_IDS = info.Patient;

Crit_Base = {};
Crit_Anes = {};
Frequ = {};
ID = {};

%% loop over all particiopants and stepsizes and calculate dpli
for p = 1:length(P_IDS)
    p_id = P_IDS{p};
    
    fprintf("Analyzing complexity of '%s' in '%s' \n", p_id,'Base');
    
    part_dir = strcat('data_aperiodic/sub-',p_id,'/eeg/');
    participant_in_Base = strcat('sub-',p_id,'_task-Base_eeg.set');
    participant_in_Anes = strcat('sub-',p_id,'_task-Anes_eeg.set');

    %% Load data
    EEG_Base = pop_loadset('filename',participant_in_Base, ...
        'filepath',part_dir);

    EEG_Anes = pop_loadset('filename',participant_in_Anes, ...
        'filepath',part_dir);

    % Resample the data
    EEG_Base = pop_resample( EEG_Base, 250);
    EEG_Anes = pop_resample( EEG_Anes, 250);
    fs = EEG_Anes.srate;

    % set up python env to use foof in the next function
    % This is for some reason not compatible with Anaconda! 
    % Run this line with your path to python if pyenv does not have a
    % version
    % pe = pyenv('Version','/Users/charlotte/Documents/Python_env/default_env/bin/python');     
    % add pytho folder to path
    % you need to pip install fooof anf numpy

    pyenv
    py.importlib.import_module('numpy')
    py.importlib.import_module('fooof')

    p_ID_Base_k = {};
    p_ID_Base_fr = {};
    p_ID_Anes_k = {};
    p_ID_Anes_fr = {};

    %set parameters
    f_range = [1,6];
    max_freq = 6;

   
    %% calculate Complexity Baseline
    parfor ch =1:size(EEG_Base.data,1)
        % find lowpass frequency
        [fr,ctr] = select_low_pass_freq(EEG_Base.data(ch,:),fs,f_range,max_freq)
        if ~isnan(fr)
            % Filter the data at fr lowpass
            EEG_Base_filt = pop_eegfiltnew(EEG_Base, [], fr, [], false, [], 0); % Lowpass filter
            % do Chaos test
            k = chaos_test(EEG_Base_filt.data(ch,:),'minmax')
            disp('Done ' + string(p_id) + ' Baseline channel ' + string(ch))
        else
            k = NaN; 
        end      
        p_ID_Base_k = [p_ID_Base_k, k];
        p_ID_Base_fr = [p_ID_Base_fr, fr];
    end

    %% calculate the same for Anesthesia
    parfor ch =1:size(EEG_Anes.data,1)
        % find lowpass frequency
        [fr,ctr] = select_low_pass_freq(EEG_Anes.data(ch,:),fs,f_range,max_freq)
        if ~isnan(fr)
            % Filter the data at fr lowpass
            EEG_Anes_filt = pop_eegfiltnew(EEG_Anes, [], fr, [], false, [], 0); % Lowpass filter
            % do Chaos test
            k = chaos_test(EEG_Anes_filt.data(ch,:),'minmax')
            disp('Done ' + string(p_id) + ' Anesthesia channel ' + string(ch))
        else
            k = NaN; 
        end      
        p_ID_Anes_k = [p_ID_Anes_k, k];
        p_ID_Anes_fr = [p_ID_Anes_fr, fr];
    end

    T_Base = array2table([cell2mat(p_ID_Base_k)', ... 
                  cell2mat(p_ID_Base_fr)'], ...
                  'VariableNames', { 'k_Base','fr_Base'});

    T_Anes = array2table([cell2mat(p_ID_Anes_k)', ... 
                  cell2mat(p_ID_Anes_fr)'], ...
                  'VariableNames', { 'k_Anes','fr_Anes'});

    % Write data to text file
    writetable(T_Base, 'Criticality_AllChannels_Base_' + string(p_id) + '.txt', "Delimiter",",")
    writetable(T_Anes, 'Criticality_AllChannels_Anes_' + string(p_id) + '.txt', "Delimiter",",")

    Crit_Base = [Crit_Base, median(cell2mat(p_ID_Base_k))];
    Crit_Anes = [Crit_Anes, median(cell2mat(p_ID_Anes_k))];
    ID = [ID, p_id];


    %% save data
    T = table(ID(:), Crit_Base(:), Crit_Anes(:), ...
        'VariableNames', { 'ID', 'Crit_Base','Crit_Anes'});
    % Write data to text file
    writetable(T, 'FINAL_Criticlity.txt')
end


T = table(ID(:), Crit_Base(:), Crit_Anes(:), ...
    'VariableNames', { 'ID', 'Crit_Base','Crit_Anes'});
% Write data to text file
writetable(T, 'Joint_Criticlity.txt')
