function k = chaos_test(y,ds_method)

% This code estimates the chaoticity of a time-series y on a scale from
% zero to one, where zero indicates periodicity/stability and one indicates
% chaos/instability.See Toker, Sommer, and D'Esposito, "A Simple Method for
% Detecting Chaos In Nature" (2020), Comm. Biol., for more details on
% time-discretization, signal normalization, and the modification to the
% 0-1 chaos test
%
% Note that for our paper, neural signals were low-pass filtered before 
% applying this analysis (see select_low_pass_freq.m)

% If y is over-sampled from a continuous-time system, take its local minima
% and maxima to transform it into a discrete-time signal
if (max(y)-min(y) )/mean(abs(diff(y))) >10
    if strcmp(ds_method,'minmax')
        y=minmaxsig(y);
    elseif strcmp(ds_method,'ds')
        oversample_flag=1;
        while oversample_flag==1
            y=downsample(y,2);
            if (max(y)-min(y) )/mean(abs(diff(y))) < 10 ...
                    || length(y)<100
                oversample_flag=0;
            end
        end
    end
end

% Normalize the standard deviation of y to 0.5
norm_fac=.5/std(y);
y=y.*norm_fac;
%y=y(1:70);
% Calculate the K-statistic using the modified 0-1 chaos test
k = z1test(y,0.5);
end



function [kmedian,p,q]=z1test(x,sig)

% The following is based on code originally written by Paul Matthews
% and made available here:
% https://www.mathworks.com/matlabcentral/fileexchange/25050-0-1-test-for-chaos
%
% This code implements one modification to Paul Matthews's code.
%
% The modification is based on Dawes and Freeland (2008), which
% includes adding a noise term to M(n). The noise is scaled by the variable
% "sig". This modification term improves distinguishability between chaotic
% vs quasi-periodic and strange non-chaotic systems

if nargin<2
    sig=1;
end

s=size(x);
if s(2)==1
    x=x';
end
N=length(x); 
j=[1:N];
t=[1:round(N/10)];
M=zeros(1,round(N/10));
c=pi/5+rand(1,1000)*3*pi/5;

for its=1:1000
    p=cumsum(x.*cos(j*c(its)));
    q=cumsum(x.*sin(j*c(its)));
    for n=1:round(N/10)
        M(n)=mean((p(n+1:N)-p(1:N-n)).^2 + (q(n+1:N)-q(1:N-n)).^2 )- ...
            mean(x)^2*(1-cos(n*c(its)))/(1-cos(c(its)))+sig*(rand-.5);
    end
    kcorr(its)=corr(t',M');
end
kmedian=median(kcorr);
end

function vec = minmaxsig(a)

inds1=islocalmax(a);
inds2=islocalmin(a);
inds=logical(inds1+inds2);
vec=a(inds);
end