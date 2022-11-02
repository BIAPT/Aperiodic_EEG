function c = criticality(k,alpha)

% This function estimates the proxmity of a system to the edge-of-chaos
% critical point, based on k (the output of chaos_test.m) and the
% parameter alpha. See Eq. 15 of the paper

if alpha<=0 || alpha>=1
    error('alpha must lie between zero and one');
end

% The k-statistic of the 0-1 chaos test can sometimes produce values
% slightly below zero or above 1. In these cases, set to either zero or one 
k(k<0)=0;
k(k>1)=1;

% Calculate proximity to edge-of-chaos criticality
c=zeros(size(k));
c(k<alpha) = k(k<alpha)./alpha;
c(k>=alpha) = 1-(k(k>=alpha)-alpha)./(1-alpha);