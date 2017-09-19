function [D, dt, t0, DataType] = anamedian(DS, Chan, DT, iCond, iRep, Twin);
% Dataset/anamean - mean over reps of analog data
%    D = anamean(DS, ChanDT, DT) returns the median of AD data from channel 
%    Chan. The median is taken after detrending individual reps using a
%    window of DT ms. No Detrending is done when DT==inf.
%    The k-th column of D contains the madiean waveform 
%    over the repetitions of the k-th stimulus condition (see stimPresent). 
%    Stimulus  conditions for which no recording is present (due to 
%    interruption of  the recording) yield Nan-valued columns.
%    See Dataset/anadata for the channel specification Chan.
%
%    D = anamedian(DS, Chan, DT, iCond) only selects the stimulus conditions
%    with indices in array iCond. This call is formally equivalent to an
%    unrestricted Anamean call followed by D=D(iCond), but the restricted
%    call is usually faster, because the nonrequested data are not read
%    from disk nor averaged. Specifying iCond values outside the range of 
%    stimulus conditions defined for DS is an error, but valid conditions 
%    for which no data are avalaible (due to interruption of the recording)
%    result in Nan-valued columns.
%
%    D = anamedian(DS, Chan, DT, iCond, iRep) also restricts the data to
%    the
%    repetitions with indices iRep. Specify iCond=0 to include all
%    stimulus conditions, iRep=0 to include all reps.
%
%    D = anamedian(DS, Chan, DT, iCond, iRep, [t0 t1]) further restricts the
%    data to the time segments starting at t0 ms and ending at t1. These 
%    boundaries must lie inside the interval [0 Dur], where Dur is the 
%    duration of the recording of a single repetition. Specifing [] is the 
%    same as not specifying the time window, resulting in the entire 
%    duration to be returned.
%
%    [D, dt, t0, DataType] = anamedian(...) also returns the sample period dt
%    in ms; time offset t0 (as determined by calibration) in ms; and
%    a char string DataType identifying the type of data stored in D.
%
%    See also Dataset/anamean, StimPresent, Dataset/anadata, Dataset/anavar,
%    Dataset/ana_myfun.

CacheFile = [mfilename '_' name(DS.ID.Experiment)];
[iCond, iRep, Twin] = arginDefaults('iCond/iRep/Twin', 0, 0, []);
iRep = replaceMatch(iRep,{[], 0, ':'}, 0); % standardize "all reps" value to avoid duplicate caches
CacheParam = {DS.ID.iDataset, DT, Chan, iCond, iRep, Twin};
Y = getcache(CacheFile, CacheParam);
if isempty(Y), % compute it and store in cache
    % real work is delegated to ana_myfun
    % === ana_myfun(DS, @(d)median(quickHP(dt,d,5),2), P.ADchan, P.iCond, P.iRep, P.Twin);
    dt = samperiod(anachan(DS, Chan));
    [D, dt, t0, DataType] = ana_myfun(DS, @(d)median(quickHP(dt,d,DT),2), Chan, iCond, iRep, Twin);
    Y = CollectInStruct(D, dt, t0, DataType);
    putcache(CacheFile, 100, CacheParam, Y);
else % unpack
    [D, dt, t0, DataType] = deal(Y.D, Y.dt, Y.t0, Y.DataType);
end







