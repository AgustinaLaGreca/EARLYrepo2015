function [pval cvs, z] = rayleigh_test(alpha)
% rayleigh_test - Rayleigh uniformness test for circular data
% [pval, cvs, z] = rayleigh_test(alpha)
%   Computes Rayleigh test for non-uniformity of circular data.
%   H0: the population is uniformly distributed around the circle
%   HA: the populatoin is not distributed uniformly around the circle
%   Assumption: the distribution has maximally one mode and the data is 
%   sampled from a von Mises distribution!
%
%   Input
%     alpha  sample of angles in cycles. If a mtrix, Each column is treated
%            as a separate dataset.
%
%   Outputs
%     pval  p-value of Rayleigh's test
%     cvs   complex vector strength
%     z     value of the z-statistic
%
% Adapted from Philippe Berens' circ_rtest (version 7/6/2008). Uses cycles
% instead of radians. Also accepts matrix data whose columns are data sets.
% This speeds up AnaRayleigh. 
% Removed binned weights and bin size arguments.
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


cvs =  mean(exp(2*pi*i*alpha)); % mean vector in complex plane
r = abs(cvs); % vector strength
n = size(alpha,1); % # samples per columns (= dataset)
% compute Rayleigh's R (equ. 27.1)
R = n*r;
% compute Rayleigh's z (equ. 27.2)
z = R.^2 / n;
% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4*n+4*(n^2-R.^2))-(1+2*n));

% outdated version:
% compute the p value using an approximation from Fisher, p. 70
% pval = exp(-z);
% if n < 50
%   pval = pval * (1 + (2*z - z^2) / (4*n) - ...
%    (24*z - 132*z^2 + 76*z^3 - 9*z^4) / (288*n^2));
% end
