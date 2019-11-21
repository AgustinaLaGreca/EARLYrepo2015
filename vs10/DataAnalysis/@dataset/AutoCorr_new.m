function [Template, GUI] = AutoCorr_new(D, figh, P)
%AutoCorr  calculates Shifted Auto Corrilogram
%   T = AutoCorr(ds) 
%   E.g.:
%           ds = dataset('A0242', '8-2');
%           T = AutoCorr(ds);
% 
%Ramses de Norre
% 
%% ---------------- CHANGELOG ----------------
% 
%  Jul 2019 Gowtham         - Adapted EvalSAC (!only for BPN!), uses just the
%                             autocorrelogram on time domain;
%                             Doesn't work realtime for now
% 
%  Thu Jan 13 2011  Abel    - Increased spkrate by factor 1000
%							- Updated userdata try/catch block
%  Mon Jan 17 2011  Abel
%	- Output in single vectors instead of repetitions. Old format is still
%	given as second output option (for now).
%
%  Tue Jan 18 2011  Abel
%   - Added option to skip THR info gathering (for fake datasets)
%
%  Fri Jan 28 2011  Abel
%   - Added sorted indepval column to output struct
%
%  Tue Feb 8 2011  Abel
%   - Added NOT normalised sac (nonorm.max(n))
%	- Added "coincidence rate (CR)" = Same as DriesNorm (CI) but not
%	normalised by square of the rate. (ac.cr(n))
%	Name changes in output:
%		ac.max -> ac.ci
%		ds.spkrate -> ac.spkrate
%		ds.nspike -> ac.nspike
%%

%% -----------------------------------template-----------------------------
Template.ds.filename     = '';        %Datafile name for datasets
Template.ds.icell        = NaN;       %Cell number of datasets
Template.ds.iseq         = NaN;       %Sequence number of first dataset
Template.ds.seqid        = '';        %Identifier of first dataset
Template.ds.isubseq      = NaN;       %Subsequence number of spiketrain used for first dataset

%Template.ds.nspike      = NaN;
%Template.ds.spkrate      = NaN;
Template.ac.nspike      = NaN;
Template.ac.spkrate      = NaN;

Template.tag             = 0;         %General purpose tag field
Template.createdby       = mfilename; %Name of MATLAB function that generated the data
Template.stim.burstdur   = NaN;       %Stimulus duration in ms
Template.stim.repdur     = NaN;       %Repetition duration in ms
Template.stim.nrep       = NaN;       %Number of repetitions
Template.stim.spl        = NaN;       %Sound pressure level in dB
Template.stim.LowFreq    = NaN;       %LowFreq for RCN
Template.stim.HighFreq   = NaN;       %HighFreq for RCN

% Template.ac.max          = NaN;
Template.ac.ci		 = NaN;	      %Maximum of shuffled autocorrelogram (DriesNorm: Correlation Index)
Template.ac.cr		 = NaN;		  %Maximum of shuffled autocorrelogram (Coincidence Rate)
Template.ac.nonorm.max   = NaN;		  %Maximum of shuffled autocorrelogram (No normalisation)

Template.ac.saczero      = NaN;       %Value at delay zero of shuffled autocorrelogram (DriesNorm)
Template.ac.peakratio    = NaN;       %Ratio of secundary versus primary peak in autocorrelogram
Template.ac.hhw          = NaN;       %Half height width on autocorrelogram (ms)
Template.ac.fft.df       = NaN;       %Dominant frequency in autocorrelogram (Hz)
Template.ac.fft.bw       = NaN;       %Bandwidth (Hz)
Template.thr.cf          = NaN;       %Characteristic frequency retrieved from threshold curve
Template.thr.sr          = NaN;       %Spontaneous rate retrieved from threshold curve
Template.thr.thr         = NaN;       %Threshold at characteristic frequency
Template.thr.q10         = NaN;       %Q10 retrieved from threshold curve
Template.thr.bw          = NaN;       %Width 10dB above threshold (Hz)
Template.rcn.thr         = NaN;       %Noise threshold (dB) derived from NSPL curve ...

%% -------------------------------default parameters-----------------------
%Calculation parameters ...
DefParam.anwin         = [0 +Inf];   %in ms (Inf designates stimulus duration)
DefParam.corbinwidth   = 0.03;       %in ms ...
DefParam.cormaxlag     = 15;         %in ms ...
DefParam.acfftrunav    = 5000;       %in Hz ...
DefParam.acfftcutoff   = 10000;       %in Hz
DefParam.calcdf        = NaN;        %in Hz, NaN (automatic), 'cf' or 'df' ...
%Plot parameters ...
DefParam.plot          = 'yes';      %'yes' or 'no' ...
DefParam.thrplot       = 'no';      %'yes' or 'no' ...
DefParam.dotplot       = 'no';      %'yes' or 'no' ...
DefParam.fftyunit      = 'dB';       %'dB' or 'P' ...
DefParam.ismashed      = false;      % Whether the dataset is mashed.
DefParam.ignorethr     = true;      % Don't look up any THR info. This option is needed for fake datasets (example: output of MergeDS)
% Yi-Hsuan
DefParam.Scarmble      = false;      %for rising floor of correlation between ACF and pusle train

% recursive call using default plot params
if numel(D) > 1 
    clf;
    for ii=1:numel(D),
        AutoCorr_new(D(ii));
    end
    return;
end

%===========single D from here=============

% handle the special case of parameter queries. Do this immediately to 
% avoid endless recursion with dataviewparam.
if isvoid(D) && isequal('params', figh),
    [H,G] = local_ParamGUI;
    return;
end

% open a new figure or use existing one?
if nargin<2 || isempty(figh),
    open_new = isempty(get(0,'CurrentFigure'));
    figh=gcf; 
else,
    open_new = isSingleHandle(figh);
end

% parameters
if nargin<3, P = []; end
if isempty(P), % use default paremeter set for this dataviewer
    P = dataviewparam(mfilename); 
end

% delegate the real work to local fcn
local_AutoCorr(D, figh, open_new, P);

% enable parameter editing when viewing offline
if isSingleHandle(figh, 'figure'), enableparamedit(D, P, figh); end;

end

function local_AutoCorr(D, figh, open_new, P)

% prepare plot
if isSingleHandle(figh, 'figure')
    figure(figh); clf; ah = gca;
    if open_new, placefig(figh, mfilename, D.Stim.GUIname); end % restore previous size 
else
    ah = axes('parent', figh);
end

Pres = D.Stim.Presentation;
P = struct(P); P = P.Param;
isortPlot = P.iCond(P.iCond<=Pres.Ncond); % limit to actual Ncond
if isortPlot==0, isortPlot = 1:Pres.Ncond; end;
Ncond = numel(isortPlot);
AW = P.Anwin;

Spt=ds.Spiketimes;

Chan = 1; % digital channel
Nrep = NrepRec(D);
SPT = spiketimes(D, Chan, 'no-unwarp');
BurstDur = GenericStimparams(D,'BurstDur');

for icond=1:Ncond, 
    if isequal('burstdur', AW)
        bdur = max(BurstDur(icond,:)); % burst dur in ms
        aw = [0 bdur]; 
    else
        aw = AW;
    end
    spt = AnWin(SPT(icond, :),aw);
    Nsp(icond) = numel([spt{:}]); 
    Rate(icond) = 1e3*Nsp(icond)./Nrep(icond)./bdur;
end; 


end

function [T,G] = local_ParamGUI % Taken from revcor
% Returns the GUI for specifying the analysis parameters.
P = GUIpanel('AutoCorr','');
iCond = ParamQuery('iCond', 'iCond:', '0', '', 'integer',...
    'Condition indices for which to calculate the CV. 0 means: all conditions.', 20);
Anwin = ParamQuery('Anwin', 'analysis window:', 'burstdur', '', 'anwin',...
    'Analysis window (in ms) [t0 t1] re the stimulus onset. The string "burstdur" means [0 t], in which t is the burst duration of the stimulus.');
maxDelay = ParamQuery('maxDelay', 'max delay:', '20', '', 'integer',...
    'Maximum delay (in ms) considered in computing the reverse correlation. 0 means: t, in which t is the burst duration of the stimulus.',1);
P = add(P, iCond);
P = add(P, Anwin, below(iCond));
P = add(P, maxDelay, below(Anwin));
P = marginalize(P,[4 4]);
G = GUIpiece([mfilename '_parameters'],[],[0 0],[10 10]);
G = add(G,P);
G = marginalize(G,[10 10]);
% list all parameters in a struct
T = VoidStruct('iCond/Anwin/maxDelay');
end