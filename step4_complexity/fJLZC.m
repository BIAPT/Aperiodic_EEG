function [norm_LZC, LZC] = fJLZC(data,fs,multivar_type, norm_type)

% This function calculates the normalized Lempel-Ziv complexity of either a
% univariate or multivariate signal. The default method for a multivariate
% signal is Zozor et al's joint Lempel-Ziv complexity (from Zozor et al,
% "On Lempel-Ziv complexity for multidimensional data analysis," 2005,
% Physica A). But, if you set multivar_type to 'concat', this function will
% instead concatenate a multivariate signal into a single univariate signal
% (following the method described in Schartner et al, "Complexity of
% Multi-Dimensional Spontaneous EEG Decreases during Propofol Induced
% General Anaesthesia," 2015, PLoS One). The LZC can either be normalized
% using a random shuffle of the binary sequence (Schartner) or a phase-
% randomized surrogate (Toker, 2022). 

% Matlab code by Daniel Toker and Diego M. Mateos. C++ code by Steeve Zozor

% Inputs
%   data - a univariate or multivariate time-series, where each row is a
%   channel and each column is a time-point 
%
%   fs - sampling frequency
%
%   multivar_type - set to 'concat' to concatenate a multivariate signal
%   into a univariate one
%
%   norm_type - type of normalizeation can be 'shuffle' or 'phase-rand'
%
% Outputs
%   norm_LZC - univariate Lempel-Ziv complexity, normalized using phase-
%   randomized Fourier transform surrogates (following method in paper)
%
%   LZC - non-normalized Lempel-Ziv complexity

% If c++ code hasn't been compiled, compile
%if exist('complexity.mexmaci64','file')~=1
%    mex complexity.c
%end

% Generate Fourier transform surrogates
if nargin==1 || isempty(fs)
    fs=1;
end

if strcmp(norm_type, 'phase-rand') 
    for i = 1:size(data,1)
        rand_data(i,:) = surrogate(data(i,:), 1, 'FT', 0, fs);
    end
end

% Binarize the input signal by thresholding at the mean of its
% instantaneous amplitude
%data=detrend(data')'; % detrend
%data=bsxfun(@minus, data,mean(data,2)); % demean
data=abs(hilbert(data'))'; % calculate instantaneous amplitude
bin_data=bsxfun(@gt,data,mean(data,2)); % binarize

if strcmp(norm_type, 'phase-rand') 
    % Binarize the phase-randomized signal by thresholding at the mean of its
    % instantaneous amplitude
    %rand_data=detrend(rand_data')'; % detrend
    %rand_data=bsxfun(@minus, rand_data,mean(rand_data,2)); % demean
    rand_data=abs(hilbert(rand_data'))'; % calculate instantaneous amplitude
    bin_rand_data=bsxfun(@gt,rand_data,mean(rand_data,2)); % binarize
end

% Shuffle not implemented for zorzor test
if size(data,1)>1 && strcmp(norm_type, 'shuffle') && ~strcmp(multivar_type, 'concat')
    disp("Method Shuffle not implemented for zozor LZC. Please use Concatenated or Joint LZC instead")
    return
end

% If using Schartner et al's concatenation method, collapse multivariate
% signal into a univariate signal
if strcmp(multivar_type, 'concat')
    bin_data=bin_data(:)';
end

% If useing shuffle normalization shuffle order of binary time series
if strcmp(norm_type, 'shuffle') 
    % get a shuffeled version of the binarized data
    bin_rand_data = bin_data(randperm(length(bin_data)));
end

% If using Schartner et al's concatenation method, collapse random 
% multivariate signal into a univariate signal
if nargin>2 && strcmp(multivar_type, 'concat')
    bin_data=bin_data(:)';
    bin_rand_data=bin_rand_data(:)';
end

%%% Joint Discretization
[Xjd, Nc]=JoinDisc(bin_data);
[rand_Xjd, rand_Nc]=JoinDisc(bin_rand_data);

%%% Lempel Ziv complexity
NAlphabet=2^Nc;
rand_NAlphabet=2^rand_Nc;
LZC=complexity(Xjd,NAlphabet);
rand_LZC=complexity(rand_Xjd,rand_NAlphabet);
norm_LZC = LZC/rand_LZC;
end

function [Xjd, Nc]=JoinDisc(M)
% Generate the joint discretization implemented by Zozor in "On Lempel-Ziv
% complexity for multidimensional data analysis"
% M is a matrix consisting of channels x time-points
[Nc, Ns]=size(M);
Vexp=(0:Nc-1);
V2=(ones(1,Nc)*2).^Vexp;
Xjd=zeros(1,Nc);


for itime=1:Ns    
    Xv=M(:,itime);
    Xjd(itime)=sum(Xv.*V2');
end
end