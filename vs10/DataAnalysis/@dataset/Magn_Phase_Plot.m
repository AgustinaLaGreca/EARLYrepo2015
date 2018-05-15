function [Y, Z]=Magn_Phase_Plot(DS, ADchan, Twin, Cdelay, varargin);
% dataset/Magn_Phase_Plot - plot magnitude & phase plots of analog data
%    S = Magn_Phase_plot(D, ADchan, Twin, Cdelay) plots magnitude and phase 
%    of analog tone responses in analog channel ADchan of D and returns the
%    result of the analysis in struct S.
%    ADchan is either a full field name of the DS.Data struct, or a number
%    1|2, which serve as shorthands for RX6_analog_1 and RX6_analog_2.
%    Twin is the analysis window, i.e., the time interval [t_start t_end] 
%    in ms over which  magnitude and phase are evaluated. The char string
%    'BurstDur' is shorthand for [0 Bd], where Bd is the burst duration of
%    the individual stimulus conditions. This is the default value. 
%    The phases are corrected for the DA/AD delay of the anaolg recording
%    (see adc_data/timelag). An additional delay correction is specified by
%    the input argument Cdelay. Default value is Cdelay = 0 ms.
%
%    [P, H] = Magn_Phase_Plot(dataset(), 1, '-paramstruct') returns a parameter
%    struct P containing all possible input parameters and their default
%    values. H is a struct containing the same fieldnames as P. Its fields
%    are brief help texts that explain the meaning of the corresponding
%    parameters.
%
%    P = Magn_Phase_Plot(DS, Chan, P), where P is a parameter struct as
%    returned by the type of call described in the previous item, uses the
%    parameters specified in P. Unspecified parameters (absent fields) are
%    set to their default values.
%
%    P = Magn_Phase_Plot(DS, Chan, 'Foo1', Value1, 'Foo2', Value2, ...)
%    is an alternative way of specifying the parameters. Parameter names,
%    here indicated as Foo1 ..., may be abbreviated as long they are
%    unambiguous and are case-insensitive (see keywordMatch).
%
%    See also dataset/powerspec.

[Twin, Cdelay] = arginDefaults('Twin/Cdelay', 'BurstDur', 0);
[Pfac, H] = local_defaultParam(ADchan); % factory defaults of params in struct format

if ischar(Twin), % if Twin resembles 'BurstDur', turn it to 'BurstDur' exactly. Same for '-paramstruct'
    TwinBD = keywordMatch(Twin, {'BurstDur' '-paramstruct'});
    if ~isempty(TwinBD), Twin = TwinBD; end
end

if isequal('-paramstruct', Twin), % return paramstruct and helper
    [Y, Z] = local_defaultParam(ADchan);
    return;
end

Pfac.doPlot = nargout<1; % unless explicitly specified, plot only when no output args are requested.
if isstruct(Twin),
    [Param, OFN] = combineStruct(Pfac, Twin, 'amend');
    if ~isempty(OFN),
        error(['Invalid parameters: ' char(10) '  ' cell2words(OFN, [char(10) '  '] )]);
    end
elseif isequal('BurstDur', Twin) || isnumeric(Twin), 
    % combine Twin with Cdelay and then with varargin, which is a series of prop/val pairs. Convert to struct
    Param = local_struct(Pfac, 'Twin', Twin, 'Cdelay', Cdelay, varargin{:});
% elseif ischar(Twin),
%     Param = local_struct(Pfac, 'Twin', 'BurstDur', 'Cdelay', 0, varargin{:});
else % (Twin,Cdelay) are in fact the 1st prop/val pair. Merge them with varargin.
    varargin = [Twin Cdelay varargin];
    Param = local_struct(Pfac, varargin{:});
end

%Param

%if Param.doPlot, clf; end
if numel(DS)>1 && has2varparams(DS(1)), % multiple RFs
    Y = [];
    for ii=1:numel(DS),
        y = Magn_Phase_Plot(DS(ii), ADchan, Twin, Cdelay, varargin{:});
        Y = [Y, y];
    end
    return;
elseif numel(DS)==1 && has2varparams(DS), % special case: single DS with two varied params
    Nx = DS.Stim.Ncond_XY(1); 
    Ny = DS.Stim.Ncond_XY(2);
    for iy=1:Ny,
        p = Param;
        p.iCond = (1:Nx) + Nx*(iy-1);
        % delegate the analysis to local fcn
        Y(iy) = local_analyze(DS, p, H);
        if Param.doPlot, 
            [ah1 ah2 hll(iy)] = local_plot(gcf, DS, Y(iy), H, lico(iy), ploma(iy));
        end
        Par2 = DS.Stim.Presentation.Y;
        LegStr{iy} = sprintf(Par2.FormatString, Par2.PlotVal(p.iCond(1)));
    end
else, % straightforward case
    for ii=1:numel(DS),
        % delegate the analysis to local fcn
        Y(ii) = local_analyze(DS(ii), Param, H);
        if Param.doPlot,
            [ah1 ah2, hll(ii)] = local_plot(gcf, DS(ii), Y(ii), H, lico(ii), ploma(ii));
        end
        LegStr{ii} = ['Rec ' num2str(DS(ii).ID.iDataset)];
    end
end
% finish plot niceties
if Param.doPlot,
    if isequal(0, Y(1).Param.Cdelay);
        phasedelayslider(ah2, 'kHz');
    else,
        hl = legend(['Cdelay = ' num2str(Param.Cdelay) ' ms'], 'location', 'NorthEast');
        lPos = get(hl,'posi'); set(hl,'pos', lPos.*[1 1.15 1 1]);
    end
    axes(ah1);
    legend(hll, LegStr, 'location', 'northwest');
    MultSets = numel(DS)>1 || isequal('RF', stimtype(DS));
    if MultSets, xlog125([0.05 50]); end
end



% 

%====================================================================
function [P, H] = local_defaultParam(ADchan);
% default parameter set with factory defaults and corresponding help texts
P.ADchan = ADchan;
H.ADchan = 'Analog channel to be analyzed. Either 1|2 or full fieldname of DS.Data.';
P.Twin = 'BurstDur';
H.Twin = 'Analysis window. Either [t0 t1] in ms or one of the strings ''BurstDur'', ''RepDur''.';
P.Cdelay = 0;
H.Cdelay = 'Compensatory delay for phase analysis.';
P.iRep = 0;
H.iRep = 'Selection of reps. 0 means all reps in dataset.';
P.iCond = 0;
H.iCond = 'Selection of stim conditions. 0 means all conditions in dataset.';
P.FreqName = 'Fcar';
H.FreqName = 'Field of DS.Stim containing the frequency array for the analysis. Default is Fcar.';
P.FreqDAC = 1;
H.FreqDAC = 'DA channel from which to extract the frequency.';
P.AvMode = 'median';
H.AvMode = 'Mode of combining the reps (char string). Either Mean or Median. Default is Median.';
P.TauDetrend = inf;
H.TauDetrend = 'Detrending time constant in ms. Default is inf ms, i.e., no detrending.';
H.FreqDAC = 'DA channel from which to extract the frequency.';
P.iHarm = 1;
H.iHarm = 'Harmonic number of the frequency to be analyzed. Default: 1, i.e., fundamental.';
P.DoNormalize = false;
H.DoNormalize = 'Whether to normalize the gain re the SPL of the stimulus component.';
P.GainShift = 0;
H.GainShift = 'Overall shift to be added to gain.';
P.doPlot = true;
H.doPlot = 'To plot or not to plot.';
P.reference = transfer();
H.reference = 'Reference against which to evaluate magn&phase. Transfer object. Default: void transfer = none.';
if ~isequal(fieldnames(P), fieldnames(H)),
    error('P and H must have same field names.');
end

function Param = local_struct(Pfac, varargin);
% convert prop/val pairs into struct after checkeing validity of props
if rem(numel(varargin),2), 
    error('Property/Values must come in pairs.');
end
Param = Pfac;
Props = varargin(1:2:end);
if ~iscellstr(Props),
    error('Parameter names (properties of prop/val pairs) must be char strings.');
end
AllProps = fieldnames(Pfac);
for ii=1:numel(Props),
    [prop, Mess] = keywordMatch(Props{ii}, AllProps, 'analysis parameter');
    error(Mess);
    Param.(prop) = varargin{2*ii};
end

function S = local_analyze(DS, P, H);
[S, CFN, CP] = getcache([mfilename '_' name(DS.Stim.Experiment)], {DS.ID P});
if ~isempty(S), return; end
% if ~ismember(stimtype(DS), {'FS', 'RF'}),
%     error(['Magn_phase_plot cannot process ' stimtype(DS) ' recordings.']);
% end
Dataset_ID = DS.ID;
ExpName = name(DS.ID.Experiment);
iDataset = DS.ID.iDataset;
DSdescr = IDstring(DS);
StimulusType = stimtype(DS);
% read data; averaged over reps
P.Twin = replaceMatch(P.Twin, 'BurstDur', []); % [] is anamean's name for burstdur 
switch lower(P.AvMode),
    case 'mean',
        [Y, dt] = anamean(DS, P.ADchan, num2cell(P.iCond), P.iRep, P.Twin);
        Y = quickHP(dt, Y, P.TauDetrend);
    case 'median',
        [Y, dt] = anamedian(DS, P.ADchan, P.TauDetrend, num2cell(P.iCond), P.iRep, P.Twin);
end
[Scal, Yunit] = conversionfactor(anachan(DS,P.ADchan));
Freq = DS.Stim.(P.FreqName);
icond = replaceMatch(P.iCond,0,':');
Freq = Freq(icond, P.FreqDAC);
IndepVarX = DS.Stim.Presentation.X;
IndepVarY = DS.Stim.Presentation.Y;
Xval = IndepVarX.PlotVal(icond);
Freq = SameSize(Freq, Xval)*P.iHarm;
if has2varparams(DS),
    fldnm = IndepVarY.FieldName;
    Yval = DS.Stim.(fldnm)(icond);
else,
    Yval = [];
end
% fit cosine
for ii=1:numel(Freq),
    Z(ii,1) = cosfit(dt, Scal*Y{ii}, Freq(ii));
    RayleighAlpha(ii,1) = anarayleigh(local_truncate(Y{ii},dt,Freq(ii)), dt, Freq(ii), 3e3/Freq(ii), 10);
end
if ~isvoid(P.reference),
    Z = Z./eval(P.reference, Freq,1);
end
Magn_dB = A2dB(abs(Z)) + P.GainShift;
if isequal('Pa', Yunit),
    Magn_dB = Magn_dB + 94;
    Yunit_dB = 'dB SPL';
elseif isequal('V', Yunit),
    Yunit_dB = 'dBV';
else,
    Yunit_dB = ['dB re 1 ' Yunit];
end
if P.DoNormalize && isequal('RF', StimulusType),
    Magn_dB = Magn_dB - DS.Stim.SPL(icond);
    Yunit_dB = 'dB';
elseif P.DoNormalize && isequal('FS', StimulusType),
    Magn_dB = Magn_dB - unique(DS.Stim.SPL);
    Yunit_dB = 'dB';
end
AD_DA_lag = timelag(anachan(DS, P.ADchan));
Phase_cycle = cangle(Z) - DS.Stim.WavePhase;
Phase_cycle = delayPhase(Phase_cycle, 1e-3*Freq, AD_DA_lag +P.Cdelay, 2);
Param = P;
%==============should be replaced by truncation in time domain============
% Remove onset delay
OnsDel = GenericStimparams(DS,'OnsetDelay'); % all onset delays
if P.iCond ~= 0
    OnsDel = OnsDel(P.iCond); % onset delay for current conditions
end
for ii = 1:numel(OnsDel)
    Phase_cycle(ii) = delayPhase(Phase_cycle(ii), 1e-3*Freq(ii), OnsDel(ii), 2);
end
Phase_cycle = cunwrap(Phase_cycle);
%=========================================================================
if isequal(1,Param.iHarm),
    Xlabel = [Param.FreqName ' (Hz)'];
else,
    Xlabel = [num2str(Param.iHarm) 'x' Param.FreqName ' (Hz)'];
end
S = CollectInStruct(ExpName, iDataset, Dataset_ID, DSdescr, StimulusType, Param, IndepVarX, IndepVarY, ...
    '-', Freq, Xval, Yval, Magn_dB, Phase_cycle, RayleighAlpha, Xlabel, Yunit_dB, '-', H);
putcache(CFN, 200, CP, S);

function [ah1 ah2 hl] = local_plot(figh, DS, Y, H, ColorSpec, PlotMarker);
figure(figh);
set(figh,'units', 'normalized', 'position', [0.28 0.38 0.62 0.52]);
ah1 = subplot(2,1,1);
PLM = pmask(Y.RayleighAlpha<=0.01);
xplot(1e-3*Y.Xval, Y.Magn_dB, '.:', ColorSpec);
hl = xplot(1e-3*Y.Xval, Y.Magn_dB+PLM, '-', 'marker', PlotMarker, ColorSpec, 'markersize', 5);
title([IDstring(DS, 'full') '  AD-'  num2str(Y.Param.ADchan) ' (' dataType(anachan(DS,Y.Param.ADchan)) ')'],...
     'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');
ylabel(['magnitude (' Y.Yunit_dB ')'],'fontsize',10);
xlog125;
set(ah1, 'XLimSpec', 'Tight'); % x-axis tight %Jan 2018
ah2 = subplot(2,1,2);
xplot(1e-3*Y.Xval, Y.Phase_cycle, '.:', ColorSpec);
xplot(1e-3*Y.Xval, Y.Phase_cycle+PLM, '-', 'marker', PlotMarker, ColorSpec);
xlabel('frequency (kHz)','fontsize',10);
set(ah2, 'XLimSpec', 'Tight'); % x-axis tight %Jan 2018
% linkaxes([ah1 ah2], 'x');
ylabel('phase (cycle)','fontsize',10);


function Y=local_truncate(Y, dt, Freq);
% truncate waveform Y so that its length matches as closely as possible an
% integer # of cycles of frequency Freq (Hz). Sample period dt in ms.
Dur = numel(Y)*dt; % ms total duration of Y
T = 1e3/Freq; % ms period of Freq
TruncDur = T*floor(Dur/T); % truncated duration. Contains integer # cycles.
NsamTrunc =round(TruncDur/dt);
Y = Y(1:NsamTrunc);







