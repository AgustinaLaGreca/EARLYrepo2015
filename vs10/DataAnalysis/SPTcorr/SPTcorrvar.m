function [V, Tau, M] = SPTcorrvar(spt1, spt2, varargin)
% SPTcorrvar - variance of spike time correlogram
%    [V, Tau] = SPTCORR(SPT1, SPT2, ...), where SPT1 and SPT2 
%    are cell arrays of spiketime arrays, returns a bootstrap estimate of
%    the variance of the correlograms produced by SPTcorr. V is the
%    variance and Tau are the bin centers of the correlogram. 
%    Ellipses ... denote the remaining parameters of the correllogram (see SPTcorr).
%
%    Each cell of SPT1 and SPT2 is considered as the response to a single 
%    "stimulus repetition". The variance is estimated by visting each
%    possible pair of reps from SPT1 and SPT2.
%
%    SPTCORR(SPT1, 'nodiag', ...) evaluates the variance across
%    non-identical pairs of the cells of SPT1.
%
%    [V, Tau, X] = SPTCORR(SPT1, ...), also returns the correlogram X
%    itself; this saves you a separate call to SPTcorr.
%
%
%    See also SPTcorr, ANWIN, XCORR.

% special recursive cases

N1 = numel(spt1);
if iscell(spt2),
    N2 = numel(spt2);
    XC = [];
    for i1=1:N1,
        for i2=1:N2,
            [xc, Tau] = SPTCORR([spt1{i1}], [spt2{i2}], varargin{:});
            XC = [XC; xc(:).']; % rows are single-pairs correlograms
        end
    end
elseif isequal('nodiag', spt2),
    XC = [];
    for i1=1:N1,
        for i2=1:N1,
            if i1~=i2,
                [xc, Tau] = SPTCORR([spt1{i1}], [spt1{i2}], varargin{:});
                XC = [XC; xc(:).']; % rows are single-pairs correlograms
            end
        end
    end
else,
    error('Both spt1 and spt2 most be cell arrys unless spt2 equals the char string ''nodiag.''');
end
V = var(XC,1); % variance across pairs
M = mean(XC,1); % mean across pairs = grand correlogram






