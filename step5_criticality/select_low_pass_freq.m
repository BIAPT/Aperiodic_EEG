function [fr,ctr] = select_low_pass_freq(x,fs,f_range,max_freq)

% This code uses the FOOOF algorithm ("Fitting Oscillations and One Over
% F") to find channel-specific slow (1-6 Hz) oscillation frequencies. To
% use, download the Matlab wrapper for the FOOOF toolbox:

% https://github.com/fooof-tools/fooof_mat


% Input
% x - a neural time-series recording
% fs - sampling frequency of the recording

% Output
% fr - the high-frequency end of the bandwidth of the slowest identified 
%      oscillation. In our paper, signals were low-pass filtered at fr

settings=struct();
[psd, freqs] = pwelch(x, [], [], [], fs);
freqs = freqs';
psd = psd';

fooof_results = fooof(freqs, psd, f_range, settings);
if ~isempty(fooof_results.gaussian_params)
    peaks = fooof_results.gaussian_params(:,1);
    inds=find(peaks<max_freq);
    if ~isempty(inds)
    for i = 1:length(inds)
        all_fr(i) = fooof_results.gaussian_params(inds(i),1)+.5*fooof_results.gaussian_params(inds(i),3);
        filt_x=eegfilt(x,fs,0,all_fr(i));
        sig_corr(i) = corr2(x,filt_x);
    end
    fr_ind = find(sig_corr==max(sig_corr)); % best match to original signal
    
    fr = fooof_results.gaussian_params(fr_ind,1)+.5*fooof_results.gaussian_params(fr_ind,3);
    ctr = fooof_results.gaussian_params(1,1);
    else
        fr=NaN;
    end
else
    fr = NaN; % return NaN if no oscillations were identified in the 1-6 Hz range
    ctr = NaN;
end