function [ds, Spt, Info, StimParam, Param] = EvalSacParseArgs(DefParam, varargin)

% Checking input arguments ...

if (length(varargin) < 1)
    error('Wrong number of input arguments.');
end
if ~isa(varargin{1}, 'dataset')
    error('First argument should be dataset.');
end

ds = varargin{1};
ParamIdx = 2;

%% Retrieving properties and checking their values ...
Param = checkproplist(DefParam, varargin{ParamIdx:end});
CheckParam(Param);

%% Unraveling and reorganizing input arguments ...
iSubSeq = 1:ds.Ncond;
%IndepVal = ds.IndepVar.Values;

Spt = ds.spiketimes(iSubSeq, :);
if DefParam.Scarmble, Spt = ScrambleSpkTr(Spt); end
ds = EmptyDataset(ds);

%% Assembling dataset information ...

Info.ds.filename = ds.expname;
Info.ds.icell    = ds.ID.iCell;
Info.ds.iseq     = ds.ID.iDataset;
Info.ds.seqid    = ds.IDstring;
Info.ds.Exp      = ds.ID.Experiment;
Info.ds.mashed   = '';
if Param.ismashed
%TODO   Info.ds.mashed = ['mashed#' num2str(ds.ID.mashed)];
end

% Collecting and reorganizing stimulus parameters ...
for n = 1:ds.Ncond
    StimParam(n) = GetStimParam(ds, iSubSeq(n));
end

%% Substitution of shortcuts in properties ...
if isinf(Param.anwin(2))
    Param.anwin(2) = min([StimParam.burstdur]);
end

%% Format parameter information ...
if isnan(Param.calcdf)
    CalcDFStr = 'auto';
elseif ischar(Param.calcdf)
    CalcDFStr = lower(Param.calcdf);
else
    CalcDFStr = EvalSacParam2Str(Param.calcdf, 'Hz', 0);
end
s = sprintf('AnWin = %s', EvalSacParam2Str(Param.anwin, 'ms', 0));
s = strvcat(s, sprintf('BinWidth = %s', EvalSacParam2Str(Param.corbinwidth, 'ms', 2)));
s = strvcat(s, sprintf('MaxLag = %s', EvalSacParam2Str(Param.cormaxlag, 'ms', 0)));
s = strvcat(s, sprintf('Calc. DF = %s', CalcDFStr));
Param.str = s;
end

%% CheckParam
function CheckParam(Param)

%Calculation parameters ...
if ~isnumeric(Param.anwin) | (size(Param.anwin) ~= [1,2]) ...
        | ~isinrange(Param.anwin, [0, +Inf]) %#ok
    error('Invalid value for property anwin.');
end
if ~isnumeric(Param.corbinwidth) || (length(Param.corbinwidth) ~= 1) ...
        || (Param.corbinwidth <= 0)
    error('Invalid value for property corbinwidth.');
end
if ~isnumeric(Param.cormaxlag) || (length(Param.cormaxlag) ~= 1) ...
        || (Param.cormaxlag <= 0)
    error('Invalid value for property cormaxlag.');
end
if ~isnumeric(Param.acfftrunav) || (length(Param.acfftrunav) ~= 1) ...
        || (Param.acfftrunav < 0)
    error('Invalid value for property acfftrunav.');
end
if ~(isnumeric(Param.calcdf) && ((Param.calcdf > 0) || isnan(Param.calcdf)))...
        && ~(ischar(Param.calcdf) && any(strcmpi(Param.calcdf, {'cf', 'df'}))), 
    error('Property calcdf must be positive integer, NaN, ''cf'' or ''df''.'); 
end

%Plot parameters ...
if ~any(strcmpi(Param.plot, {'yes', 'no'}))
    error('Property plot must be ''yes'' or ''no''.');
end
if ~(Param.acfftcutoff > 0)
    error('Invalid value for property acfftcutoff.');
end
if ~any(strcmpi(Param.fftyunit, {'dB', 'P'}))
    error('Property fftyunit must be ''dB'' or ''P''.');
end
end
%% StimParam
function StimParam = GetStimParam(ds, iSubSeq)

%If stimulus and repetition duration is not the same for different channels or different 
%datasets or both, then the minimum value is used. NaN's are automatically removed by min
%function ...
StimParam.burstdur = round(burstdur(ds,iSubSeq));
StimParam.repdur   = round(ds.repdur(:));

%Repetition numbers for both datasets ...
StimParam.nrep = ds.NrepRec(iSubSeq);

%Sound pressure level (dB) is always returned as a two element vector ... When one of the
%datasets is a binaural dataset and the SPL administered at the ears is different, then 
%the SPLs are combined using the power average ...
spl1 = GetSPL(ds);
StimParam.spl = CombineSPLs(spl1(iSubSeq, 1), spl1(iSubSeq, end));

StimParam.LowFreq = ds.stimparam.LowFreq;
StimParam.HighFreq = ds.stimparam.HighFreq;

% Gowtham, July 2019: right now this is just for BPN, could be adapted for other functions
if(isfield(ds.Stim,'CutoffFreq'))
    StimParam.CutoffFreq=ds.Stim.CutoffFreq(iSubSeq);
end

%Format stimulus parameters ...
s = sprintf('BurstDur = %s', EvalSacParam2Str(StimParam.burstdur, 'ms', 0));
s = strvcat(s, sprintf('IntDur = %s', EvalSacParam2Str(StimParam.repdur, 'ms', 0)));
s = strvcat(s, sprintf('#Reps = %s', EvalSacParam2Str(StimParam.nrep, '', 0)));
s = strvcat(s, sprintf('SPL = %s', EvalSacParam2Str(unique(StimParam.spl), 'dB', 0)));
s = strvcat(s, sprintf('LowFreq = %s', EvalSacParam2Str(unique(StimParam.LowFreq), 'Hz', 0)));
s = strvcat(s, sprintf('HighFreq = %s', EvalSacParam2Str(unique(StimParam.HighFreq), 'Hz', 0)));
StimParam.str = s;
end

%--------------------------------------------------------------------------------------------------
function SpkOut = ScrambleSpkTr(SpkIn)
% shuffle spike train by randomizing the ISI's
if iscell(SpkIn)
    SpkOut = cell(size(SpkIn));
    for i = 1:size(SpkIn,1)
        for j = 1:size(SpkIn,2)
            SpkOut{i,j}= ScrambleSpkTr(SpkIn{i,j});
        end
    end
else
    NInt = length(SpkIn);
    D = diff([0 sort(SpkIn)]); % first sort to make sure diff gives the ISI's
    SpkOut = cumsum(D(randperm(NInt)));
end
end
%--------------------------------------------------------------------------------------------------
