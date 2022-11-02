%% Charlotte Maschke June 2022
% This script goal is to calculate complexity measures with the code published by Toker et.al 
% https://figshare.com/articles/software/Consciousness_is_supported_by_near-critical_cortical_electrodynamics/12949355

% Source Setup
%
INPUT_DIR = 'Users/charlotte/Documents/GitHub/one_over_f/data_aperiodic/';
OUTPUT_DIR = 'RESULTS/criticality/';

% load participant info
opts = detectImportOptions("data_aperiodic/data_2states.txt",'Delimiter','\t');
info = readtable("data_aperiodic/data_2states.txt",opts);
P_IDS = info.Patient;

Chaos_Base = {};
Crit_Base = {};
Freq_Base = {};
Nopeak_Base = {};

Chaos_Anes = {};
Crit_Anes = {};
Freq_Anes = {};
Nopeak_Anes = {};
ID = {};

%% loop over all particiopants and stepsizes and calculate Criticality
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
    p_ID_Base_nopeak = 0;
    total_Base = 0;
    p_ID_Base_fr = {};
    p_ID_Anes_k = {};
    p_ID_Anes_nopeak = 0;
    total_Anes = 0;
    p_ID_Anes_fr = {};

    %set parameters
    f_range = [1,6];
    max_freq = 6;

    data_Base = EEG_Base.data;
    data_Anes = EEG_Anes.data;

    % split data into non-overlaping 10 s windows
    window_size = 10; % in seconds
    step_size = 10; % in seconds
    [Base_epochs] = create_window(data_Base, fs, window_size, step_size);
    [Anes_epochs] = create_window(data_Anes, fs, window_size, step_size);
   
    % take same nr of trials, maximal 30
    nr_trials = min([size(Base_epochs,1),size(Anes_epochs,1),30]);

    %% calculate Complexity Baseline
    for trial = 1:nr_trials
        data_trial = squeeze(Base_epochs(trial,:,:));
        tmp_ch_Base = {}; % collect data for trial
        parfor ch =1:size(EEG_Base.data,1)
            % find lowpass frequency
            data_channel = data_trial(ch,:);
            [fr,ctr] = select_low_pass_freq(data_channel,fs,f_range,max_freq)
            if ~isnan(fr)
                % Filter the data at fr lowpass
                data_channel_filt=eegfilt(data_channel,fs,0,fr);
                % do Chaos test
                k = chaos_test(data_channel_filt,'minmax')
                p_ID_Base_fr = [p_ID_Base_fr, fr];
            else
                k = NaN; 
                p_ID_Base_nopeak = p_ID_Base_nopeak + 1
            end    
            tmp_ch_Base = [tmp_ch_Base, k];
            total_Base = total_Base + 1
        end
        p_ID_Base_k = [p_ID_Base_k, median(cell2mat(tmp_ch_Base),'omitnan')];
        disp('Done ' + string(p_id) + ' Baseline Trial ' + string(trial))
    end

    %% calculate the same for Anesthesia
    for trial = 1:nr_trials
        data_trial = squeeze(Anes_epochs(trial,:,:));
        tmp_ch_Anes = {}; % collect data for trial
        parfor ch =1:size(EEG_Anes.data,1)
            % find lowpass frequency
            data_channel = data_trial(ch,:);
            [fr,ctr] = select_low_pass_freq(data_channel,fs,f_range,max_freq)
            if ~isnan(fr)
                % Filter the data at fr lowpass
                data_channel_filt=eegfilt(data_channel,fs,0,fr);
                % do Chaos test
                k = chaos_test(data_channel_filt,'minmax')
                p_ID_Anes_fr = [p_ID_Anes_fr, fr];
            else
                k = NaN; 
                p_ID_Anes_nopeak = p_ID_Anes_nopeak + 1
            end    
            tmp_ch_Anes = [tmp_ch_Anes, k];
            total_Anes = total_Anes + 1
        end
        p_ID_Anes_k = [p_ID_Anes_k, median(cell2mat(tmp_ch_Anes),'omitnan')];
        disp('Done ' + string(p_id) + ' Anesthesia Trial ' + string(trial))
    end

    alpha =  0.85; % as defined by Toker et al 2022
    k_Base = median(cell2mat(p_ID_Base_k));
    c_Base = criticality(k_Base,alpha);

    k_Anes = median(cell2mat(p_ID_Anes_k));
    c_Anes = criticality(k_Anes,alpha);

    Chaos_Base = [Chaos_Base, k_Base ];
    Freq_Base = [Freq_Base, median(cell2mat(p_ID_Base_fr))]; % mean over all frequencies
    Crit_Base = [Crit_Base, c_Base];
    % percentage of nopeak from all channels
    Nopeak_Base = [Nopeak_Base, (p_ID_Base_nopeak/total_Base)*100 ];

    Chaos_Anes = [Chaos_Anes, k_Anes];
    Crit_Anes = [Crit_Anes, c_Anes];
    Freq_Anes = [Freq_Anes, median(cell2mat(p_ID_Anes_fr))]; % mean over all frequencies
     % percentage of nopeak from all channels
    Nopeak_Anes = [Nopeak_Anes, (p_ID_Anes_nopeak/total_Anes)*100 ];

    ID = [ID, p_id];

    %% save data
    T = table(ID(:), Chaos_Base(:), Chaos_Anes(:), Crit_Base(:), Crit_Anes(:), Freq_Base(:), Freq_Anes(:), Nopeak_Base(:), Nopeak_Anes(:),...
        'VariableNames', { 'ID', 'Chaos_Base','Chaos_Anes', 'Crit_Base','Crit_Anes','Freq_Base','Freq_Anes','Nopeak_Base','Nopeak_Anes'});
    % Write data to text file
    writetable(T, 'EOC_Criticality.txt')
end


T = table(ID(:), Chaos_Base(:), Chaos_Anes(:), Crit_Base(:), Crit_Anes(:),Freq_Base(:), Freq_Anes(:), Nopeak_Base(:), Nopeak_Anes(:),...
    'VariableNames', { 'ID', 'Chaos_Base','Chaos_Anes', 'Crit_Base','Crit_Anes','Freq_Base','Freq_Anes','Nopeak_Base','Nopeak_Anes'});
% Write data to text file
writetable(T, 'EOC_Criticality.txt')
