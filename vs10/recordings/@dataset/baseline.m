function [D, dt, t0, DataType] = baseline(DS, Chan, Flag, Twin);
% Dataset/baseline - extract baseline analog data from Dataset
%    D = anadata(DS, Chan, 'pre') returns the AD data from channel Chan
%    that correspond to the baseline interval preceding the stimulus. 
%
%    D = baseline(DS, Chan, 'post') returns the AD data from channel Chan
%    that correspond to the baseline interval trailing the stimulus.
%
%    baseline(DS, Chan, Flag, [t0 t1]) only returns a restricted portion 
%    ranging from t0 to t1 ms, where t=0 refers to the start of the 
%    requested baseline interval.
%
%    [D, dt, t0, DataType] = baseline(...) also returns the sample period 
%    dt in ms; time offset t0 (as determined by calibration) in ms; and
%    a char string DataType identifying the type of data stored in D.
%
%    See also StimPresent, Dataset/anadata, adc_data/samples.

[Twin] = arginDefaults('Twin', []); % default: max time window

AD = anachan(DS, Chan); % appropriate field of DS.Data
StimPres = DS.Stim.Presentation;


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

switch lower(Flag),
    case 'pre',
        ioffset = 0;
        nmax = StimPres.NsamPre;
    case 'post',
        ioffset = StimPres.SamOffset(end-1);
        nmax = StimPres.NsamPost;
    otherwise,
        error('Invalid flag. Must be ''pre'' or ''post''.');
end
if isequal(0,nmax),
end

if isnan(iAnwin1), % no time window; get whole rep
    NsamRead = nmax;
else, % check & apply time window
    if iAnwin1>nmax,
        error('Time window Twin exceeds baseline duration.');
    end
    NsamRead = iAnwin1-iAnwin0;
    ioffset = ioffset + iAnwin0;
end
[D, dt, t0, DataType] = samples(AD, ioffset+1, ioffset+NsamRead);










