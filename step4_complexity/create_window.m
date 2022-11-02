% This function is to get overlapping windowed data
function [windowed_data] = create_window(data, sfreq, window, step)
    % input data needs to be n_channels x time
    % step: stepsize in seconds
    % window size in seconds
    % sfreq sampling frequency

    length_recording = size(data,2);
    number_channels = size(data,1);

    window_size = window * sfreq; % in points
    step_size = step*sfreq;
    iterator = 1:step_size:(length_recording - window_size);
    windowed_data = zeros(length(iterator),number_channels,window_size);
    index = 1;
    for i = 1:step_size:(length_recording - window_size)
        windowed_data(index,:,:) = data(:,i:i+window_size-1);
        index = index + 1;
    end
end