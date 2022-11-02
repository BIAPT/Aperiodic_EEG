function stochastic = stochastic_test(surr_y,fs)

% Time-series test of signal stochasticity. The test has been modified from
% its original version, as described in Toker, Sommer, and D'Esposito, "A
% simple method for detecting chaos in nature," Comms. Bio. (2020).

% The test calculates the permutation entropy of a time-series and
% compares it to the permutation entropies of 1,000 Amplitude Adjusted
% Fourier Transform surrogates and 1,000 Cyclic Phase Permuation
% surrogates. If the permutation entropy of the original signal falls
% within either surrogate distribution, the signal is classified as
% stochastic. The test has been modified from its original form to catch a
% failure case of the Cyclic Phase Permutation surrogate algorithm, which
% is that it can generate surrogates with identical permutation entropies
% for some low-noise or noise-free signals. If this failure case occurs,
% "jitter" in the form of 2.5% white noise is incrementally added to the 
% original signal until this failure case is broken. 

% Input
% surr_y - a recorded time-series signal
% fs - the sampling frequency of the recorded signal

% Output
% stochastic - will equal 1 for a predominantly stochastic signal, and will
%              equal 0 for a predominantly deterministic signal

if size(surr_y,1)>size(surr_y,2)
    surr_y=surr_y';
end

% Generate 1,000 Amplitude Adjusted Fourier Transform surrogates and
% calculate their permutation entropies
try
    [surr, params] = surrogate(surr_y, 1000, 'AAFT', 1, fs);
catch
    [surr, params] = surrogate(zscore(surr_y), 1000, 'AAFT', 1, fs);
end
sig=params.cutsig;
for i = 1:1000
    surr_h1(i) = petropy(surr(i,:),8,1);
end
perm_h1 = petropy(sig,8,1);

stochastic1 = perm_h1>=min(surr_h1)&&perm_h1<=max(surr_h1);

if stochastic1==1
    stochastic=1;
else
    
    % Stochastic *nonlinear* data might pass the prior test and be classified
    % as deterministic. To try to rule this out, do a second test of determinism using
    % cyclic phase permutation surrogates, which maintains individual cycles in
    % the data but shuffles them, thus breaking determinism
    surr_flag=0;
    n_level = 0.025;
    orig_surr_y=surr_y;
    while surr_flag==0
        try
            [surr, params] = surrogate(surr_y, 1000, 'CPP', 1, fs);
        catch
            [surr, params] = surrogate(zscore(surr_y), 1000, 'CPP', 1, fs);
        end
        sig=params.cutsig;
        for i = 1:1000
            x1 = surr(i,:);
            surr_h2(i) = petropy(x1,8,1);
        end
        % If surrugoates don't have unique permutation entropies, add 2.5%
        % noise
        if length(unique(surr_h2))<500
            surr_y=orig_surr_y+randn(1,length(surr_y)).*n_level.*std(orig_surr_y,[],2);
            n_level=n_level+0.025;
        else
            surr_flag=1;
        end
    end
     
    perm_h2 = petropy(sig,8,1);
    stochastic = double(perm_h2>=min(surr_h2)&&perm_h2<=max(surr_h2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code for calculating permutation entropy is written by
% Andreas Muller, and is available on the University of Potsdam's TOCSY -
% Toolboxes for Complex Systems: http://tocsy.pik-potsdam.de/petropy.php

function H = petropy(x,n,tau,method,accu)
% The petropy function calculates the permutation entropy of data series.
%
% version 1.1, 27.02.2015:
%   - corrected typo in the description ('same' -> 'equal')
%   - line 106: changed unique-function in newer MATLAB-version
%
% Permutation Entropy
%
% H = PETROPY(X,N,TAU,METHOD,ACCU) computes the permutation entropy H of
% a scalar vector X, using permutation order N, time lags from TAU and
% METHOD to treat equal values. The ACCU parameter describes the accuracy
% of the values in X by the number of decimal places.
%
% x      - data vector (Mx1 or 1xM)
% n      - permutation order
% tau    - time lag scalar OR time lag vector (length = n-1)
% method - method how to treat equal values
%   'noise' - add small noise
%   'equal' - allow same rank for equal values
%   'order' - consider order of appearance (first occurence --> lower rank)
% accu   - maximum number of decimal places in x
%         (only used for method 'noise')
%
% References:
%
% Bandt, C.; Pompe, B. Permutation Entropy: A Natural Complexity
% Measure for  Time Series. Phys. Rev. Lett. 88 (2002) 17, 174102
%
% Riedl, M.; Müller, A.; Wessel, N.: Practical considerations of
% permutation entropy. The European Physical Journal Special Topics
% 222 (2013) 2, 249?262
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example:
%
% H = petropy([6,9,11,12,8,13,5],3,1,'order');
% H =
%       1.5219
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5, accu = 4; end
if nargin < 4, method = 'order'; end

x = x(:);
M = length(x);
equal = false;

if n*log10(n)>15, error('permutation dimension too high'); end
if (length(tau)) > 1 && (length(tau) ~= n-1), error('time lag vector has to have n-1 entries'); end
if ((n-1)*min(tau) >=M) || max(tau) >= M, error('too few data points for desired dimension and lags'); end


switch lower(method)
    case {'noise'}
        %disp('Method: add small noise')
        x = x + rand(M,1)*10^(-accu-1);
    case 'equal'
        %disp('Method: allow equal ranks')
        equal = true;
    case 'order'
        %disp('Method: consider order of occurrence')
    otherwise
        error('unknown method')
end

if length(tau) > 1
    tau = reshape(tau,length(tau),1);
    tau = sort([0;tau]);
    % build n x (M-tau(n))-matrix from shifted values in x
    shift_mat = zeros(n,M-tau(n));
    for ii=1:n
        shift_mat(ii,:) = x(tau(ii)+1:M-tau(n)+tau(ii));
    end
else
    % vectorized
    shift_mat_ind = reshape(0:tau:(n-1)*tau,[],1) * ones(1,M-(n-1)*tau) +...
        ones(n, 1) * reshape(1:(M-(n-1)*tau),1,[]);
    shift_mat = x(shift_mat_ind);
end

if equal
    % allow equal values the same index
    ind_mat = zeros(size(shift_mat));
    for ii=1:size(ind_mat,2)
        [~,~,ind_mat(:,ii)]=unique(shift_mat(:,ii),'first');
    end
else
    % sort matrix along rows to build rank orders, equal values retain
    % order of appearance
    [~, sort_ind_mat] = sort(shift_mat,1);
    ind_mat = zeros(size(sort_ind_mat));
    for ii=1:size(ind_mat,2)
        ind_mat(sort_ind_mat(:,ii),ii) = 1:n;
    end
end
% assign unique number to each pattern (base-n number system)
ind_vec = n.^(0:n-1) * (ind_mat-1);

% find first occurence of unique values in 'ind_vec' and use
% difference to determine length of sequence of the same numbers; e.g.
% sort_ind_vec = [21 21 11 19 11], unique_values = [11 19 21],
% ia = [1 3 4]: 11 occurs on places #1 and #2, 19 on #3 and 21 on #4 and #5
[~,ia,~] = unique(sort(ind_vec), 'first');

% use the following line, if you are using MATLAB in a lower version than
% 8.3:
% permpat_num = diff([ia (length(ind_vec)+1)]);
permpat_num = diff([ia; (length(ind_vec)+1)]);

permpat_num = permpat_num/sum(permpat_num);

% compute permutation entropy
H = -sum(permpat_num .* log2(permpat_num));
