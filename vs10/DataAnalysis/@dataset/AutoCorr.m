function [Template, GUI] = AutoCorr(D, figh, P)
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
        AutoCorr(D(ii));
    end
    return;
end

%===========single D from here=============

% handle the special case of parameter queries. Do this immediately to 
% avoid endless recursion with dataviewparam.
if isvoid(D) && isequal('params', figh),
    [Template, GUI] = local_ParamGUI;
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
% if isempty(P), % use default paremeter set for this dataviewer
%     P = dataviewparam(mfilename); 
% end

Ncond = D.Stim.Presentation.Ncond;
Nrep = D.Stim.Presentation.Nrep;

%% main program

% Evaluate input arguments ...
if (nargin == 1) && ischar(D) && strcmpi(D, 'factory')
    disp('Properties and their factory defaults:')
    disp(DefParam);
    return;
else
    [ds, Spt, Info, StimParam, Param] = ...
        EvalSacParseArgs(DefParam, D );
end

% to display the Exp info on the last grid of the plot
PlotInfo = [Info.ds.filename " Rec: " Info.ds.iseq  ":" Info.ds.seqid ];
PlotInfo = join(PlotInfo);

% THR is not used in this particular dataviewer. [below functions need to be modified]
Rcn = struct([]);
[CF, SR, Thr, BW, Q10, Str] = deal(NaN);
Thr = lowerFields(CollectInStruct(CF, SR, Thr, BW, Q10, Str));

% Calculate data
for n = 1:Ncond
    CalcData(n) = CalcDifCor(Spt(n, :), Thr, Rcn, Info, Param, StimParam(n));
    CalcData(n).ds.isubseq = n;
    PlotDifCor(CalcData, D, Param,figh, StimParam, Ncond, PlotInfo, open_new);
end

% % Calculate data
% for m = 1:Ncond
%     for n = 2:Nrep
%         CalcData(n) = CalcDifCor(Spt(n, :), Thr, Rcn, Info, Param, StimParam(n));
%         CalcData(n).ds.isubseq = n;
%         PlotDifCor(CalcData, D, Param,figh, StimParam, Ncond, PlotInfo, open_new);
%     end
% end

% PlotDifCor(CalcData, ds, Param,figh, StimParam);

end

%% CalcData
function CalcData = CalcDifCor(Spt, Thr, Rcn, Info, Param, StimParam)

WinDur = abs(diff(Param.anwin)); %Duration of analysis window in ms ...

%Correlation of noise token responses of a cell with the responses of that same
%cell to that same noise token. The spiketrains are derived from the same cell
%so this is called a Shuffled AutoCorrelogram (or SAC).
[Ysac, T, NC] = SPTCORR(anwin(Spt, Param.anwin), 'nodiag', ...
    Param.cormaxlag, Param.corbinwidth, WinDur);
spkrate = NC.Rate1;
nspike = NC.Nspike1;

%by Abel: save not normalised
CalcData.ac.nonorm.co = Ysac;
CalcData.ac.nonorm.max = max(Ysac);

%by Abel: "coincidence rate" = same as normalised to Driesnorm (correlation
%index)(see SPTCORR()) but no normalised for average rate (square rate)
NRep = StimParam.nrep;
CalcData.ac.cr = (1000*Ysac)/((NRep*(NRep-1))*WinDur);
CalcData.ac.cr = max(CalcData.ac.cr);

%Apply Driesnorm
Ysac = ApplyNorm(Ysac, NC);

%Performing spectrum analysis on the SAC. Because an autocorrelogram has a DC
%component this is removed first ...
FFTsac = spectana(T, detrend(Ysac, 'constant'), 'RunAvUnit', 'Hz', ...
    'RunAvRange', Param.acfftrunav);
%The magnitude spectrum of a correlogram function is actually a power spectrum,
%therefore all magnitude units need to be changed ...
FFTsac.Magn.P  = FFTsac.Magn.A;
FFTsac.Magn.A  = sqrt(FFTsac.Magn.A);
FFTsac.Magn.dB = FFTsac.Magn.dB/2;

%Determine which dominant frequency to be used in the calculation ...
if ~isempty(Thr)
    DomFreq = DetermineCalcDF(Param.calcdf, Thr.cf, NaN, FFTsac.DF);
else
    DomFreq = DetermineCalcDF(Param.calcdf, NaN, NaN, FFTsac.DF);
end

if (DomFreq ~= 0)
    DomPer = 1000/DomFreq;
else
    DomPer = NaN;
end %Dominant period in ms ...

%Calculating the ratio of the secondary over the primary peak for the SAC and
%the half height width on the peak (relative to asymptote one)...
HalfMaxSac = ((max(Ysac)-1)/2)+1;
SacHHWx = cintersect(T, Ysac, HalfMaxSac);
SacHHW = abs(diff(SacHHWx));
[SacPeakRatio, SacXPeaks, SacYPeaks] = getpeakratio(T, Ysac, DomFreq);

%Reorganizing calculated data ...
CalcData.delay         = T;
CalcData.ac.normco     = Ysac;
CalcData.ac.saczero    = Ysac(T == 0);
CalcData.ac.ci		   = max(Ysac);
CalcData.ac.peakratio  = SacPeakRatio;
CalcData.ac.peakratiox = SacXPeaks;
CalcData.ac.peakratioy = SacYPeaks;
CalcData.ac.hhw        = SacHHW;
CalcData.ac.hhwx       = SacHHWx;
CalcData.ac.halfmax    = HalfMaxSac;
CalcData.ac.fft.freq   = FFTsac.Freq;
CalcData.ac.fft.p      = FFTsac.Magn.P;
CalcData.ac.fft.db     = FFTsac.Magn.dB;
CalcData.ac.fft.df     = FFTsac.DF;
CalcData.ac.fft.bw     = FFTsac.BW;
CalcData.ds            = Info.ds;
CalcData.thr           = Thr;
%delete rcn for runing popscript
CalcData.rcn           = Rcn;
%by Abel: increase spkrate by *1000
CalcData.ac.spkrate    = spkrate * 1000;
CalcData.ac.nspike     = nspike;

% PlotDifCor(CalcData, ds, Param,figh, npans);

end

%% ApplyNorm
function Y = ApplyNorm(Y, N)
if ~all(Y == 0)
    Y = Y/N.DriesNorm;
else
    Y = ones(size(Y));
end
end

%% DetermineCalcDF
function DF = DetermineCalcDF(ParamCalcDF, ThrCF, DifDF, SacDF)

if isnumeric(ParamCalcDF)
    if ~isnan(ParamCalcDF)
        DF = ParamCalcDF;
    elseif ~isnan(DifDF)
        DF = DifDF;
    elseif ~isnan(SacDF)
        DF = SacDF;
    else
        DF = ThrCF;
    end
elseif strcmpi(ParamCalcDF, 'cf')
    DF = ThrCF;
elseif strcmpi(ParamCalcDF, 'df')
    if ~isnan(DifDF)
        DF = DifDF;
    else
        DF = SacDF;
    end
else
    DF = NaN;
end
end

%% PlotDifCor
function PlotDifCor(CalcData, ds, Param, figh, StimParam, Ncond, PlotInfo, open_new)

if isSingleHandle(figh, 'figure')
    figure(figh); clf; ah = gca;
    if open_new, placefig(figh, mfilename, ds.Stim.GUIname); end % restore previous size 
else
    ah = axes('parent', figh);
end


% Get X and Y data and convert to cell arrays
Xsac = {CalcData.delay};
Ysac = [CalcData.ac];
Ysac = {Ysac.normco};

npans=size(Xsac);
npans=npans(2);

[axh, Lh, Bh] = plotpanes(Ncond+1,1/2,figh);

% Gowtham, July 2019: right now this is just for BPN, could be adapted for other functions
if(isfield(StimParam,'CutoffFreq'))
    for i=1:npans
        CutoffFreq(i)=strcat(num2str(getfield(StimParam,{i},'CutoffFreq')), " Hz");
    end
end

% if order is by condition
for i=1:npans
    h = axh(i);
    axes(h);
    plot(h,Xsac{i},Ysac{:,i});
    if(isfield(StimParam,'CutoffFreq'))     title(h, CutoffFreq(i)); end
end
axes(axh(end)); text(0.1, 0.5, PlotInfo, 'fontsize', 12, 'fontweight', 'bold','interpreter','none');
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