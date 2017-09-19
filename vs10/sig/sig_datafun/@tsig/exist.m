function E = exist(T, interval);
% tsig/exist - time ranges at which the channels of a tsig are defined
%    E=exist(T,[t0,t1]), where T is a multi-channel tsig, returns a struct
%    array E whose elements indicate the successive time intervals within
%    the total time range (t0,t1) at which the same subset of channels of 
%    T exist (different "episodes"). From E(k) to E(k+1), some channel 
%    kicks in or disappears. The fields of E are
%         tstart: start [ms] of this interval
%           tend: end [ms] of this interval
%          ichan: index array indicating the channels that exist during the
%                 interval (tstart,tend)
%         istart: sample index corresponding to tstart, counted from t0
%           iend: sample index corresponding to tend, counted from t0
%             i0: starting sample index, counted from start of channel
%             i1: last sample index, counted from start of channel
%           nsam: # samples is current episode
%
%    E=median(T) uses t0=min(T.onset)  and t1=max(T.offset) for the
%    evaluation interval.
%    
%    EXAMPLE.
%    Let T be sampled at 10 kHz and have 2 channels, the first of which 
%    runs from 5..20 ms; the second one from 10..50 ms. Then exist(T) 
%    returns a 3x1 struct array as displayed in the following table:
%
%              tstart   tend   istart iend  existing
%       E(1)     5       10       1    50     [1]
%       E(2)    10       20       1    50     [1 2]
%       E(3)    20       50       1    50     [2]
%
%    In tsig methods that analyze things across channels, tsig/exist
%    is used in cases that do not permit default values in non-existing
%    channels. The analysis then truly depends on the subset of channels
%    existing at a given instant of time, and tsig/exist provides the
%    bookkeeping for a blocked approach. An example is tsig/median.
%
%    See also tsig/median, tsig/plus, tsig/cat, tsig/isvoid.

if ~isTsig(T); error('First argument must be tsig object.'); end
if nargin<2, interval=[min(onset(T)) max(offset(T))]; end

if ~isequal([1 2], size(interval)),
    error('Second argument (interval) must be 2-element row vector.');
end
error(numericTest(interval, '', 'Interval argument is '));

if isvoid(T),
    error('Invalid input: void tsig object.'); 
end

[t0,t1] = deal(interval(1), interval(2));

% generate single-channel tsig T0 which spans entire range (t0,t1)
DT = DT(T); % sample period in ms
t2i = @(t)round(t/DT); % portable time-to-sample converter
Nsam = t2i(t1-t0);
if Nsam<1,
    E = [];
    return;
end
T0 = tsig(Fsam(T), zeros(Nsam,1),t0);
Tonset = T.t0; % for later reference
% Expand the time range of all channels of T by adding T0 to T, but make
% sure the gaps are recognizable: fill them with NaNs, not zeros
T = binop(T,T0,@plus,nan,0);
T = [T.Waveform{:}]; % dump waveforms in big matrix
[Nsam, Nchan] = size(T);
Ina = ~isnan(T); % find out where the non-nans are ...
Ina = Ina.* repmat(2.^(0:Nchan-1),Nsam,1); % ... & in which channel they are
Ina = sum(Ina,2); % col array of bitmasks indicating existence at that time of the channels
% The "episodes" of the help text are characterized by some channel kicking
% in or out, i.e., a change in Ina.
ichange = [0; find(diff(Ina)); length(Ina)]; % the flanking 0 & length(Ina) mark the birth & death of anything
Nepi = length(ichange)-1;
for iepi=1:Nepi,
    istart = ichange(iepi)+1;
    iend = ichange(iepi+1);
    tstart = t0+DT*(istart-1);
    tend = t0+DT*(iend);
    ichan = find(~isnan(T(istart,:)));
    i0 = 1+t2i(tstart-Tonset(ichan));
    i1 = t2i(tend-Tonset(ichan));
    nsam = iend-istart+1;
    E(iepi) = CollectInStruct(tstart, tend, istart, iend, ichan,i0,i1,nsam);
end



