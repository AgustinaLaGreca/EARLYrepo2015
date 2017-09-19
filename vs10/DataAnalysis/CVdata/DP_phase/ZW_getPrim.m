function Q = ZW_getPrim(DS, iChan, iCond, ME)
%ZW_getDP -return ampl. & phase for primaries in analog zwuis data
%
% == SYNTAX ==
%   Q = ZW_getDP(DS, iChan, iCond, ME)
% == INPUTS ==
%   DS:     dataset (see getDS.m)
%   iChan:  numeric ID for the analog channel
%   iCond:  numeric ID for the stimulus condition
%   ME:     transfer object used to normalize spectra. If omitted/empty,
%           not applied. Default = [].
%
% == OUTPUT ==
%   Q:      struct with the following fields:
%
%    Fprim:    array with primary frequencies (in Hz)
%    ph:        array with phase for primaries (in cycle), 
%    mg:        array with amplitude for primaries (in dB)
%    r:         array with vector strengths for primaries 
%    wght:      array with weight factors for primaries, based on alpha (=0.001)
%               and mg > noise floor.
%    alpha:     array with conf. levels for vector strengths
%    ___________
%    Nchunk:    10: =number of chunks used to calculate r and alpha
%    rampdur:   20: =duration of ramps in calculation of r and alpha
%
%
% See also ZW_getDP, SpecComp, rayleighspec, ZWsignal, ZWspec, ZWnoise

Nchunk = 10; %see rayleighspec.m
rampdur = 20; %in ms; see rayleighspec.m

%--- handle the inputs ---
[ME] = arginDefaults('ME',[]);

%--- frequencies of primaries; in Hz ---
Fprim = DS.Stim.Fzwuis;
Fprim = Fprim(:);

%--- calculate the spectra; loopmean over 'once' the zwuis-periodicity ---
[df, Aspec, Pspec]=ZWspec(DS, iChan, iCond, ME,'once');

%--- get primaries from spectra ---
ph = SpecComp(df, Pspec, Fprim/1e3);   %phase; in cycle
mg = SpecComp(df, Aspec, Fprim/1e3);   %amplitude; in dB

%--- calc. vector strength (r) and corresponding conf. level (alpha) for all DPs ---
[dt, Y]=ZWsignal(DS, iChan, iCond, 'once');
% [r, alpha] = rayleighspec(Y, dt, Fprim, rampdur, Nchunk);
[r, alpha] = anarayleigh(Y, dt, Fprim, rampdur, Nchunk);
alpha = alpha.Alpha;

%--- find the noise floor at the primary frequencies ---
noise = ZWnoise(DS, iChan, iCond, ME, Fprim);

%--- convert noise and alpha to "digital" (0 or 1) weight factor ---
wght = zeros(size(Fprim));
iKeep = mg>noise & alpha==0.001;
wght(iKeep)=1;

%--- create output ---
Q = CollectInStruct(Fprim, ph, mg, r, wght, alpha,'-', Nchunk, rampdur);




