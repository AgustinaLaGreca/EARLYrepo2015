function Q = ZW_getDP(DS, iChan, iCond, ME, Norder, phaseref)
%ZW_getDP -return ampl. & phase for specified DP-orders in analog zwuis data
%
% == SYNTAX ==
%   Q = ZW_getDP(DS, iChan, iCond, ME, Norder, phaseRef)
% == INPUTS ==
%   DS:     dataset (see getDS.m)
%   iChan:  numeric ID for the analog channel
%   iCond:  numeric ID for the stimulus condition
%   ME:     transfer object used to normalize spectra. If omitted/empty,
%           not applied. Default = [].
%   Norder: (array of) DP-orders to include. If given as an array
%           all the specified DP-orders are returned. Default = 3.
%           NOTE: orders > 5 are not allowed.
%   phaseref: string, either 'response' or 'stimulus' which indicates
%           how the DP-phase should be referenced. 'response' refers to the
%           phase of the primary components in the recorded signal.
%           'stimulus' to the "electrical" phase of the primaries
%           (i.e. as send to the driver). Default = 'response'. Testing for
%           its value is case-insensitive and may be abbreviated.
%
% == OUTPUT ==
%   Q:      struct with the following fields:
%    Fzwuis:    array with stimulus frequencies (in Hz)
%    Fdp:       array with DP-frequencies (in Hz)
%    ph:        array with phase for DPs (in cycle), referenced according to phaseRef
%    mg:        array with amplitude for DPs (in dB)
%    r:         array with vector strengths for DPs 
%    wght:      array with weight factors for DPs, based on alpha (=0.001)
%               and mg > noise floor. Whether DP is unique or not, is not considered in wght
%    alpha:     array with conf. levels for vector strengths
%    notUnique: sparse logical array indicating whether DP-freq is unique (0) or not (1)
%    notallPrim:sparse logical array indicating whether all primaries for
%               DP were significant (0) or not (1)
%    ___________
%    Norder:    (array of) DP-orders that are considered
%    phaseRef:  string to identify phase referencing
%    Nchunk:    10: =number of chunks used to calculate r and alpha
%    rampdur:   20: =duration of ramps in calculation of r and alpha
%    ___________
%    M:         sparse matrix such that Fdp = M*DS.Stim.Fzwuis(:)
%
%
% See also DPfreqs, SpecComp, rayleighspec, ZWsignal, ZWspec, ZWnoise

Nchunk = 10; %see rayleighspec.m
rampdur = 20; %in ms; see rayleighspec.m

%--- handle the inputs ---
[ME, Norder, phaseref] = arginDefaults('ME/Norder/phaseref',[],3,'response');
[phaseref, mess] = keywordMatch(phaseref, {'response' 'stimulus'});
error(mess);

Norder = unique(Norder);
if any(Norder>5),
    error('Specified DP-orders contains value > 5. This is too much to handle...');
end

%--- frequencies of primaries; in Hz ---
Fzwuis = DS.Stim.Fzwuis;
Fzwuis = Fzwuis(:);

%--- get info on DPs ---
S = DPfreqs(Fzwuis, Norder(end),1);
S(~ismember([S.order],Norder))=[];  %restrict to specified orders (via Norder)
notUnique = [S.mult]>1;             %find the "non-unique" DPs
notUnique = sparse(notUnique(:));   %make sparse array to save space
Fdp = [S.freq]';                    %the remaining DP-freqs; in Hz

%--- create matrix M, such that Fdp = M*Fprim(:) ---
tmp = [S.weight]; tmp = tmp(:);
M = reshape(tmp,numel(S(1).weight),numel(S))';
M = sparse(M); % convert M to sparse matrix to save space.

%--- calculate the spectra; loopmean over 'once' the zwuis-periodicity ---
[df, Aspec, Pspec]=ZWspec(DS, iChan, iCond, ME,'once');

%--- prepare for application of "phaseRef" on phase data ---
if isequal(phaseref, 'response'),
    PHprim = SpecComp(df, Pspec, Fzwuis/1e3); %phase re. prims from spectrum
else,
    PHprim = DS.Stim.StartPhase; %phase re. "electrical" prims (as send to driver)
end
PHcorr = M*PHprim(:); %phase correction due to "phaseRef"

%--- get DPs from spectra ---
ph = SpecComp(df, Pspec, Fdp/1e3)-PHcorr;   %phase; in cycle
mg = SpecComp(df, Aspec, Fdp/1e3);          %amplitude; in dB

%--- calc. vector strength (r) and corresponding conf. level (alpha) for all DPs ---
[dt, Y]=ZWsignal(DS, iChan, iCond, 'once');
% [r, alpha] = rayleighspec(Y, dt, Fdp, rampdur, Nchunk);
[r, alpha] = anarayleigh(Y, dt, Fdp, rampdur, Nchunk);
alpha = alpha.Alpha;

%--- find the noise floor at the DP-frequencies ---
DPnoise = ZWnoise(DS, iChan, iCond, ME, Fdp);

%--- convert DPnoise and alpha to "digital" (0 or 1) weight factor ---
wght = zeros(size(Fdp));
iKeep = mg>DPnoise & alpha==0.001;
wght(iKeep)=1;

%--- figure out for each DP whether all primaries that went into generating
%   it are significant ---
T = DPmatrix2sum(M);
P = ZW_getPrim(DS, iChan, iCond, ME);
Pw = P.wght;
notallPrim = true(size(Fdp));

for ii = 1:numel(Fdp),
    notallPrim(ii) = any(Pw(T(ii,:))==0);
end
notallPrim = sparse(notallPrim);

%--- create output ---
Q = CollectInStruct(Fzwuis, Fdp, ph, mg, r, wght, alpha, notUnique, notallPrim,...
    '-',Norder, phaseref, Nchunk, rampdur, '-', M);




