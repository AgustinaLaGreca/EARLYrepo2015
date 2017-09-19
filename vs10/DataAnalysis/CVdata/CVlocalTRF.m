function LTR = CVlocalTRF(DSA, DSB, Chan, DX, doPlot);
% CVlocalTRF - local BM transfer fcn a la Re, He & Porsov, 2011
%    LTR = CVlocalTRF(DS1, DS2, Chan, DX, doPlot)
%    Inputs DS1 and DS2 are arrays of ZW datasets from bead A and bead B,
%    respectively. Output LTR is a struct array holding the local transfer
%    functions A->B of "matched pairs" of recordings. In order to be a
%    matching pair, DSA(I) and DSB(J) must have identical stimulus
%    frequencies, SPL and profile params (if a profile is imposed).
%    Chan is DAC channel. DX is distance between beads in um.
%
%    CVlocalTRF(LTR, 'plot'), where LTR is the output of a previous call 
%    to CVlocalTRF, plots the local transfer function(s).
%
%    CVlocalTRF(LTR, 'plotns') disables both the phaseslider and shifting 
%    of the phase curves toward zero, i.e., it uses the exact phase values
%    contained in LTR.
%
%    CVlocalTRF(LTR, 'plot', ha) and CVlocalTRF(LTR, 'plot',ha) 
%    use axes handles ha(1) and ha(2) instead of creating new ones.
%
%    CVlocalTRF(LTR, 'plotA'), where LTR is the output of a previous call 
%    to CVlocalTRF, plots the transfer of bead A.
%
%    CVlocalTRF(LTR, 'plotB'), where LTR is the output of a previous call 
%    to CVlocalTRF, plots the transfer of bead A.
%
%    CVlocalTRF
[doPlot, Chan] = arginDefaults('doPlot, Chan',1, []);

if nargin<4 && ischar(DSB),
    switch lower(DSB),
        case 'plot',
            local_plot(DSA,1, Chan);
        case 'plotns',
            local_plot(DSA,0, Chan);
        case 'plota',
            apple(dataset(), DSA.AP_A, 'plot');
        case 'plotb',
            apple(dataset(), DSA.AP_B, 'plot');
    end
    return;
end

CritStimFLDs = {'GUIname' 'LowFreq' 'HighFreq' 'Ncomp' 'BaseFreq' 'FreqTol' 'ZWseed' 'SPL' ...
    'ProfileType' 'ProfileSPLjump' 'SPLjump' 'NsamCycle'};

% extract critical stim params
[StimA, irecA] = local_stimpar(DSA, CritStimFLDs);
[StimB, irecB] = local_stimpar(DSB, CritStimFLDs);
% compare stimuli btw A and B
[iA, iB, iorphA, iorphB] = local_compareStim(StimA, StimB);
if ~isempty(iorphA),
    warning(['Orphaned recordings bead A: rec(s) ', trimspace(num2str(irecA(iorphA))) '.'] );
end
if ~isempty(iorphB),
    warning(['Orphaned recordings bead B: rec(s) ', trimspace(num2str(irecB(iorphB))) '.'] );
end
% sort stimulus par to be consistent with datasets
StimA = StimA(iA);
StimB = StimB(iB);
% apple analysis
for ii=1:numel(iA),
    dsa = DSA(iA(ii));
    dsb = DSB(iB(ii));
    if ~isequal('-', StimA(ii).ProfileType),
        error('Analysis of profiled zwuis rec not yet implemented.');
    end
    AP_A(ii) = apple(dsa, Chan, 1);
    AP_B(ii) = apple(dsb, Chan, 1);
    A2B(ii) = local_trf(AP_A(ii), AP_B(ii), DX, StimA(ii));
end
ExpName = expname(dsa);
Stim = StimA;
LTR =CollectInStruct(ExpName, irecA, irecB, Stim, AP_A, AP_B, A2B);

if doPlot,
    local_plot(LTR);
end

%==============================
function [Stim, irec] = local_stimpar(DS, CritStimFLDs);
N = numel(DS);
irec = zeros(1,N);
for ii=1:N,
    ds = DS(ii);
    if ispooled(ds),
        ds = members(ds);
        ds = ds(1);
    end
    irec(ii) = ds.irec;
    st = ds.Stim;
    Stim(ii) = structpart(st, CritStimFLDs);
end

function [iA, iB, iorphA, iorphB] = local_compareStim(StimA, StimB);
[iA, iB, iorphA, iorphB] = deal([]);
NA = numel(StimA);
NB = numel(StimB);
for ii = 1:NA,
    % find matching B
    for jj = 1:NB,
        if isequal(StimA(ii), StimB(jj));
            iA = [iA, ii];
            iB = [iB, jj];
            break;
        end
    end
end
iorphA = setdiff(1:NA, iA);
iorphB = setdiff(1:NB, iB);
stA = StimA(iA);
[iA, iB] = sortAccord(iA, iB, 1e3*[stA.LowFreq]+[stA.SPL]+0.01*[stA.ProfileSPLjump]);

function A2B = local_trf(AP_A, AP_B, DX, Stim);
ExpName = AP_A.ExpName;
irecA = AP_A.irec;
irecB = AP_B.irec;
if ~isequal(AP_A.Fprim, AP_A.Fprim)m,e
    error(sprintf('Unequal Fprim between recs %d and %d.', irecA, irecB));
end
Fprim = AP_A.Fprim;
FprimD = (AP_A.Fprim(1:end-1)+AP_A.Fprim(2:end))/2; % halfway frqs for derivative metrics like U
Gain = AP_B.Gain - AP_A.Gain;
Phase = cunwrap(AP_B.Phase - AP_A.Phase);
Alpha = max(AP_B.Alpha, AP_A.Alpha);
AlphaD = max(Alpha(1:end-1),Alpha(2:end));
PHsig = Phase(Alpha<=0.001);
if ~isempty(PHsig),
    Phase = Phase - ceil(max(PHsig));
end
c = -DX*1e-6*Fprim./Phase;
U = -DX*1e-6*diff(Fprim)./diff(Phase);
cD = interp1(Fprim, c, FprimD);
U_over_c = U./cD;
ProfileType = AP_A.ProfileType;
baseSPL = AP_A.baseSPL;
SPLjump = AP_A.SPLjump;
LegStr = AP_A.LegStr;
A2B = CollectInStruct(ExpName, irecA, irecB, DX, '-zwuis', Fprim, Gain, Phase, Alpha, ...
    '-velo', FprimD, AlphaD, c, cD, U, U_over_c, ...
    '-zwuisparams', ProfileType, baseSPL, SPLjump, LegStr);
RMS = local_RMS(AP_A, AP_B, A2B, Stim);
A2B = structJoin(A2B, RMS);

function RMS = local_RMS(APA, APB, A2B, stim);
% compute RMS of recs
PowA = dB2P(APA.Gain+stim.SPL+stim.SPLjump);
PowA(APA.Alpha>0.001) = 0;
PowB = dB2P(APB.Gain+stim.SPL+stim.SPLjump);
PowB(APB.Alpha>0.001) = 0;
RMSv_A = P2dB(sum(PowA));
RMSx_A = P2dB(sum(1e-3*APA.Fprim.*PowA));
RMSv_B = P2dB(sum(PowB));
RMSx_B = P2dB(sum(1e-3*APB.Fprim.*PowB));
RMSv_AB = 0.5*(RMSv_A+RMSv_B);
RMSx_AB = 0.5*(RMSx_A+RMSx_B);
RMS = CollectInStruct('-RMS', RMSv_A, RMSx_A, RMSv_B, RMSx_B, RMSv_AB, RMSx_AB);

%=============================================
function local_plot(LTR, doShift, ha);
[doShift, ha] = arginDefaults('doShift, ha', 1, []);
if isempty(ha),
    set(gcf,'units', 'normalized', 'position', [0.354 0.16 0.479 0.744]);
    ha(1) = subplot(2,1,1);
    ha(2) = subplot(2,1,2);
end
% ==GAIN
N = numel(LTR.A2B);
Cmap = jet(N);
axes(ha(1));
for ii=N:-1:1,
    a2b = LTR.A2B(ii);
    CLR = struct('color', Cmap(ii,:), 'markerfacecolor', Cmap(ii,:));
    xplot(a2b.Fprim/1e3, a2b.Gain,':',CLR); 
    hl(ii) = xplot(a2b.Fprim/1e3, a2b.Gain+pmask(a2b.Alpha<=0.001),'o-', CLR, 'linewidth', 2, 'markersize', 4); 
end
legend(hl, {LTR.A2B.LegStr}, 'location', 'northwest');
xlog125;
grid on;
if isfield(LTR.AP_A, 'irec_pool'),
    title([LTR(1).ExpName  char(10) numcell2str([LTR.AP_A.irec_pool]) char(10) ' vs ' numcell2str([LTR.AP_B.irec_pool])], ...
        'fontsize', 12, 'interpreter', 'none');
else,
    title(LTR(1).ExpName, 'fontsize', 14);
end
% ===PHASE
axes(ha(2));
for ii=N:-1:1,
    a2b = LTR.A2B(ii);
    CLR = struct('color', Cmap(ii,:), 'markerfacecolor', Cmap(ii,:));
    PH = a2b.Phase;
    if doShift,
        ianchor = find(a2b.Alpha<=0.001,1,'first');
        if ~isempty(ianchor),
            PH = PH - round(PH(ianchor));
        end
    end
    xplot(a2b.Fprim/1e3, PH,':',CLR); 
    hl(ii) = xplot(a2b.Fprim/1e3, PH+pmask(a2b.Alpha<=0.001),'o-', CLR, 'linewidth', 2, 'markersize', 4); 
end
xlog125;
if doShift,
    phasedelayslider(ha(2), 'kHz');
end
grid on;
linkaxes(ha, 'x');



















