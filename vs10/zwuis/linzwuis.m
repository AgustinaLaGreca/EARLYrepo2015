function S = linzwuis(Fbound, Fbase, Rseed, Nit);
% linzwuis - array of frequencies not coinciding with DPs of orders 2 or 3
%   linzwuis(Fbound, Fbase) attempts to determine a set of frequencies f(k)
%   that are all integer multiples of Fbase and obeying
%   Fbound(k,1) <= f(k) <= Fbound(k,2), for which 
%      f(k) ~= f(l) for all k,l;
%      f(k) - f(l) ~= f(m)     for all k,l,m;
%      f(k) + f(l) - f(m) ~= f(n) for all k,l,m,n with k~=l~=m~=n.
%
%   linzwuis(Fbound, Fbase, Rseed) uses random seed Rseed. The default is
%   Rseed=nan, i.e. "random random seed". See SetRandState.
%
%   linzwuis(Fbound, Fbase, Rseed, N) uses a max of N iterations to minimize
%   coinciding components.

[Rseed, Nit] = arginDefaults('Rseed/Nit', nan, 200);
% if ~isnan(Rseed), Nit = 1;
% end
Rseed = SetRandState(Rseed);

[S, CFN, CP] = getcache(mfilename, CollectInStruct(Fbound, Fbase, Rseed, Nit));
if ~isempty(S), return; end

Nf = size(Fbound,1); % # freqs
nFbound = round(Fbound/Fbase); % normalize freqs re Fbase
M = [eye(Nf); DPmatrix(Nf,2,0)]; % DPmatrix producing primaries and weightsum-zero 2nd-order DPs

MinMult = inf;
Nmult2 = inf;
for iit=1:Nit,
    % find unique set of freqs
    Nmult1 = 1;
    Nclash12 = 1;
    ism = 0;
    while Nmult1>0 || Nclash12>0 && ism<2*Nit,
        nFreq = local_randFreq(nFbound); % candidate primary freqs
        Nmult1 = numel(nFreq)-numel(unique(deciRound(nFreq,7))); % non-uniqueness of primaries
        DPs = DPfreqs(nFreq,2,0);
        Nclash12 = sum(ismember(deciRound([DPs.freq],7), deciRound(nFreq,7))); % #DP freqs coinciding w any primary
        ism = ism+1;
    end
    if Nmult1>0 || Nclash12>0, break; end
    % check for 2nd-order, sumweight-zero DPs to be unique & different from
    % all primaries
    Nmult2 = sum([DPs.mult]>1);
    if Nmult2<MinMult,
        MinMult = Nmult2;
    end
    if Nmult2==0, break; end
end

% collect return args
Nit = iit;
Freq = nFreq*Fbase; % un-normalize
    S = CollectInStruct(Fbound, Fbase, '-', Freq, Nmult2, Rseed, Nit, nFreq, '-', Nmult1, Nclash12);
putcache(CFN, 500, CP, S);

% ===========
function freq = local_randFreq(Fbound);
Nchoice = 1+Fbound*[-1; 1]; % # potential choices for each cmp
Nc = size(Fbound,1);
freq = Fbound(:,1) + floor(Nchoice.*rand(Nc,1));








