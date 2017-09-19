function [D, dt, t0, DataType] = anavar(DS, Chan, iCond, iRep, Twin);
% Dataset/anavar - variance across reps of analog data
%    D = anavar(DS, Chan) returns the mean of AD data from channel 
%    Chan. The k-th column of D contains the time-by-time variance across
%    the repetitions of the k-th stimulus condition (see stimPresent).
%    See Dataset/anadata for the channel specification Chan.
%
%    D = anavar(DS, Chan, iCond) only selects the stimulus conditions
%    with indices in array iCond. This call is formally equivalent to an
%    unrestricted anavar call followed by D=D(iCond), but the restricted
%    call is usually faster, because the nonrequested data are not read
%    from disk nor averaged.
%
%    D = anavar(DS, Chan, iCond, iRep) also restricts the data to the
%    repetitions with indices iRep. Specify iCond=0 to include all
%    stimulus conditions, iRep=0 to include all reps.
%
%    [D, dt, t0, DataType] = anavar(...) also returns the sample period dt
%    in ms; time offset t0 (as determined by calibration) in ms; and
%    a char string DataType identifying the type of data stored in D.
%
%    See also StimPresent, Dataset/anadata, Dataset/anamean, Dataset/anamean, 
%    Dataset/ana_myfun.

[iCond, iRep, Twin] = arginDefaults('iCond/iRep/Twin', 0, 0, []);
% real work is delegated to ana_myfun
[D, dt, t0, DataType] = ana_myfun(DS, @(d)var(d,0, 2), Chan, iCond, iRep, Twin);







