function [D, dt, t0, DataType] = anadata(DS, Chan, iCond, iRep, Twin);
% Dataset/anadata - extract analog data from Dataset
%    D = anadata(DS, Chan, iCond) returns the AD data from channel Chan
%    that correspond to a single stimulus condition iCond. A matrix D is 
%    returned whose columns correspond to the repetitions of the stimulus 
%    condition. If Chan is a number, say 2, data are obtained from the 
%    RX6_analog_2 field of DS.Data. Alternatively, a full field name may be 
%    specified, e.g., Chan = 'RX8_analog_13'. Stimulus conditions are 
%    denoted by their index as iCond in DS.Stim.Presentation. 
%
%    D = anadata(DS, Chan, iCond, iRep) only returns the repetitions
%    indicated by index array iRep. This is formally equivalent to an
%    unrestricted anadata call followed by D=D(:,iRep), but the restricted
%    call generally faster because the nonspecified reps are never read.
%    Specifying iRep=0 is the same as not specifying it, resulting in all
%    reps to be returned.
%
%    D = anadata(DS, Chan, iCond, iRep, [t0 t1]) further restricts the
%    data to the time segments starting at t0 ms and ending at t1. These 
%    boundaries must lie inside the interval [0 Dur], where Dur is the 
%    duration of the recording of a single repetition. Specifing [] is the 
%    same as not specifying the time window, resulting in the entire 
%    duration to be returned.
%
%    [D, dt, t0, DataType] = anadata(...) also returns the sample period dt
%    in ms; time offset t0 (as determined by calibration) in ms; and
%    a char string DataType identifying the type of data stored in D.
%
%    See also StimPresent, Dataset/anamean, Dataset/anavar,
%      adc_data/Samples, Dataset/baseline.

[iRep, Twin] = arginDefaults('iRep/Twin', 0, []); % default: all reps; max time window
iRep = replaceMatch(iRep,{[], 0, ':'}, 0); % standardize "all reps" value

if ~isscalar(iCond),
    error('iCond argument must be single condition index.');
end

AD = anachan(DS, Chan); % appropriate field of DS.Data
StimPres = DS.Stim.Presentation;

% check whether requested conditions and reps are valid
if (iCond<1) || (iCond>StimPres.Ncond);
    error('Stimulus-condition index iCond exceeds data dimensions.');
end

RepSelection = ~isequal(iRep,0); % true when specific subset of reps was specified
if RepSelection, 
    if any(iRep<1) || any(iRep>StimPres.Nrep);
        error('Stimulus-repetition index iRep exceeds data dimensions.');
    end
end
% find presentation counts ipres corresponding to given stim conditions
ipres = find(StimPres.iCond==iCond); 
nsamRep = StimPres.NsamPres(ipres(1));
% reduce ipres according to completeness of recording
ipres = ipres(ipres<=NpresRecorded(DS));
IM = ismember(StimPres.iRep(ipres), iRep);
if RepSelection, % further reduce ipres using iRep spec
    ipres = ipres(ismember(StimPres.iRep(ipres), iRep));
    if isempty(ipres),
        error('No recordings with specified iRep present in (interrupted) dataset.');
    end
else, % no specific iRep selection -> anything goes
end


D = [];
if isempty(Twin),
    [iAnwin0, iAnwin1] = deal(0,nan);
else,
    if ~isequal(2,numel(Twin)) || any(Twin<0),
        error('Invalid time window Twin.');
    end
    iAnwin0 = round(Twin(1)/samperiod(AD));
    iAnwin1 = round(Twin(2)/samperiod(AD));
end
if isempty(ipres), % return empty array + other outputs
    D = zeros(nsamRep,0);
    [dum, dt, t0, DataType] = samples(AD, 0,-1);
end
for ii=1:numel(ipres),
    nsamRep = StimPres.NsamPres(ipres(ii));
    ioffset = StimPres.SamOffset(ipres(ii));
    if isnan(iAnwin1), % no time window; get whole rep
        NsamRead = nsamRep;
    else, % check & apply time window
        if iAnwin1>nsamRep,
            error('Time window Twin exceeds Repdur');
        end
        NsamRead = iAnwin1-iAnwin0;
        ioffset = ioffset + iAnwin0;
    end
    [d, dt, t0, DataType] = samples(AD, ioffset+1, ioffset+NsamRead);
    %dsize(D,d)
    D = [D, d];
end










