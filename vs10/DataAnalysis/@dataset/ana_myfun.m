function [D, dt, t0, DataType] = ana_myfun(DS, Fun, Chan, iCond, iRep, Twin);
% Dataset/ana_myfun - apply a function to reps of analog data
%    D = ana_myfun(DS, @Foo, Chan) applies function Foo to each of the
%    matrices of AD data consisting of the reps of a stimulus condition.
%
%    D = ana_myfun(DS, @Foo, Chan, iCond) only selects the stimulus conditions
%    with indices in array iCond. This call is formally equivalent to an
%    unrestricted Ana_myfun call followed by D=D(iCond), but the restricted
%    call is usually faster, because the nonrequested data are not read
%    from disk nor averaged. 
%
%    D = ana_myfun(DS, @Foo, Chan, iCond, iRep) also restricts the data to 
%    the repetitions with indices iRep. Specify iCond=0 to include all
%    stimulus conditions.
%
%    D = ana_myfun(DS, @Foo, Chan, iCond, iRep, [t0 t1]) further restricts the
%    data to the time segments starting at t0 ms and ending at t1. These 
%    boundaries must lie inside the interval [0 Dur], where Dur is the 
%    duration of the recording of a single repetition. Specifing [] is the 
%    same as not specifying the time window, resulting in the entire 
%    duration to be returned.
%
%    [D, dt, t0, DataType] = ana_myfun(...) also returns the sample period 
%    dt in ms; time offset t0 (as determined by calibration) in ms; and
%    a char string DataType identifying the type of data stored in D.
%
%    Dataset/ana_myfun is the workhorse for Dataset/anamean, etc.
%
%    See also StimPresent, Dataset/anadata, Dataset/anamean,
%    Dataset/anavar, adc_data.

[iCond, iRep, Twin] = arginDefaults('iCond/iRep/Twin', 0, 0, []);
makeCell = iscell(iCond);
if makeCell, iCond = cell2mat(iCond); end
Pres = DS.Stim.Presentation;

iCond = replaceMatch(iCond, {0, ':'}, 1:Pres.Ncond);

if makeCell, D = {}; else, D = []; end
for ic = iCond(:).',
    [d, dt, t0, DataType] = anadata(DS, Chan, ic, iRep, Twin);
    D = [D, Fun(d)];
end







