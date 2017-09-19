function [D, dt, t0, DataType] = anastd(DS, Chan, iCond, iRep, Twin);
% Dataset/anastd - standard deviation across reps of analog data
%    D = anastd(DS, Chan) returns the standard deviation of AD data from 
%    channel Chan. The k-th column of D contains the standard deviation 
%    across corresponding time points of the repetitions of the k-th 
%    stimulus condition (see stimPresent).
%    See Dataset/anadata for the channel specification Chan.
%
%    D = anastd(DS, Chan, iCond) only selects the stimulus conditions
%    with indices in array iCond. This call is formally equivalent to an
%    unrestricted Anastd call followed by D=D(iCond), but the restricted
%    call is usually faster, because the nonrequested data are not read
%    from disk nor averaged.
%
%    D = anastd(DS, Chan, iCond, iRep) also restricts the data to the
%    repetitions with indices iRep. Specify iCond=0 to include all
%    stimulus conditions, iRep=0 to include all reps.
%
%    D = anastd(DS, Chan, iCond, iRep, [t0 t1]) further restricts the
%    data to the time segments starting at t0 ms and ending at t1. These 
%    boundaries must lie inside the interval [0 Dur], where Dur is the 
%    duration of the recording of a single repetition. Specifing [] is the 
%    same as not specifying the time window, resulting in the entire 
%    duration to be returned.
%
%    [D, dt, t0, DataType] = anastd(...) also returns the sample period dt
%    in ms; time offset t0 (as determined by calibration) in ms; and
%    a char string DataType identifying the type of data stored in D.
%
%    See also StimPresent, Dataset/anadata, Dataset/anavar,
%    Dataset/ana_myfun.

[iCond, iRep, Twin] = arginDefaults('iCond/iRep/Twin', 0, 0, []);
% real work is delegated to ana_myfun
[D, dt, t0, DataType] = ana_myfun(DS, @(d)std(d,0,2), Chan, iCond, iRep, Twin);







