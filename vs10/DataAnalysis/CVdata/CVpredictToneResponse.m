function [S,Z] = PHASE_CVpredictToneResponse(A, freq, SPL, Rcrit, tau)
% CVpredictToneResponse - predict tonal responses from zwuis responses
%    S = CVpredictToneResponse(A, freq, SPL, Rcrit, tau), with A apple
%    output from a  family of flat zwuis conditions (no middle-ear
%    corrections!), returns  a struct S in which the predictions for tonal
%    response are provided. freq may be arrays or matrices; they will be
%    "samesized" (like vector inputs of meshgrid). Rcrit is Rayleigh
%    criterium used to filter the apple out A. Only data whose Rayleigh
%    alpha<=Rcrit will enter the prediction machinery. Before
%    interpolation, phase is delayed by tau ms.
%
%    The outputS contains the fields
%
%       freq: [Hz] tone frequencies. Identical to input arg freq, but may
%              be same-sized.
%        SPL: [dB SPL] tone level Identical to input arg freq, but may
%              be same-sized.
%     SPLeqx: [dB SPL] equivalent zwuis SPL/tone according to
%             displacement-driven gain control (XGC).
%      Gainx: [dB] tone gain (response/stimulus amp) predicted
%             displacement-driven gain control (XGC).
%     Phasex: [cycle] tone phase predicted displacement-driven gain control
%       Resx: [dB] absolute mismatch in gain-controlling output power and
%             resulting output power. Metric for self-consistency of the
%             procedure.
%     SPLeqv: [dB SPL] same as SPLeqx, but based on velocity-driven gain
%             control (VGC)
%      Gainv: [dB] VGC version of Gainx
%     Phasev: [cycle] tone phase predicted velocity-driven gain control
%       Resv: [dB] VGC version of Resx
% 
%    All of these fields have the same size as determined by the call
%    [freq, SPL] = sameSize(freq, SPL).
%
%    If freq and SPL are not specified, they will be matrices covering the
%    frequencies and SPLs of the apple output A, except that the SPL will
%    range ~30 dB higher (this anticipates a downward shift of effective
%    levels).
%
%    [S, Z] = CVpredictToneResponse(A, ...) also returns struct Z
%    containing the apple output A in a readily plotted form. Phase of Z
%    has also been delayed by tau.
%
%    EXAMPLE
%       [S,Z] = CVpredictToneResponse(A);
%       subplot(2,2,1); 
%       plot(S.freq/1e3, S.Gainx); 
%       xlog125([5 30]); 
%       subplot(2,2,3);
%       plot(S.freq/1e3, S.Phasex); 
%       subplot(2,2,2); 
%       plot(Z.Fprim/1e3, Z.Gain_zw);
%       xlog125([5 30]);
%       subplot(2,2,4);
%       plot(Z.Fprim/1e3, Z.Phase_zw); 
%
%    See also dataset/apple; dataset/ampl_phase_plot.

[freq, SPL, Rcrit, tau] = arginDefaults('freq/SPL/Rcrit/tau',[],[],0.05,0);
if isempty(freq),
    freq = A(1).Fprim; 
    freq = freq(:);
end
if isempty(SPL),
    SPL = unique([A.baseSPL]);
    maxSPL = max(SPL); 
    dSPL = mean(diff(SPL));
    SPL = [SPL(:).' (maxSPL+dSPL:dSPL:maxSPL+35)];
end

[freq, SPL] = SameSize(freq, SPL);

if ~isequal('-', unique([A.ProfileType])),
    error('All apple outputs must originate from flat zwuis conditions');
end
A = sortAccord(A, [A.baseSPL]); % sort SPLs
GG = cat(1, A.Gain) + pmask(cat(1,A.Alpha)<=Rcrit);
SPLzw = cat(1, A.baseSPL);
allSPL = SameSize(SPLzw, GG);
Fprim = A(1).Fprim;
allFprim = SameSize(Fprim, GG);

% make finely spaced lookup table
fSPL = linspace(min(SPLzw),max(SPLzw),300).';
fGG = interp1(SPLzw, GG, fSPL);
ffreq = logispace(min(Fprim), max(Fprim), 200);
f2j = @(f)round(interp1(log(ffreq), 1:numel(ffreq), log(f))); % freq -> index in fGG 
fGG = interp1(Fprim.', fGG.', ffreq.').';

% Phase interpolation
PP = cat(1, A.Phase) + pmask(cat(1,A.Alpha)<=Rcrit);
PP = delayPhase(PP.', SameSize(Fprim.',PP.')./1e3, tau, 2); % delay by tau
fPP = interp1(Fprim.', PP, ffreq(1,:).').'; % linear interpolation across frequencies
cPP = exp(i*2*pi*fPP); % complex phase
fPP = cunwrap(cangle(interp1(SPLzw, cPP, fSPL))); % complex plane linear interpolation across SPLs due to plateau
% fPP = fPP-repmat(round(fPP(:,1)),1,numel(fPP(1,:)));
[fSPL, ffreq] = SameSize(fSPL, ffreq);
fPP = fPP-SameSize(round(nanmean(fPP)),fPP);

Px_out = nansum(dB2P(30+allSPL+GG)./allFprim,2); % total output power per SPL, BM displacement
Pv_out = nansum(dB2P(allSPL+GG),2); % total output power per SPL, BM velocity
fPx_dB = interp1(SPLzw, P2dB(Px_out), fSPL);
fPv_dB = interp1(SPLzw, P2dB(Pv_out), fSPL);

% visit all tonal components one by one and find the most self-consistent
% gain setting, i.e., the one whose output power gives rise to the gain
% setting that gave rise to its output power that gave rise ...
[Resx, SPLeqx, Gainx, Phasex, Resv, SPLeqv, Gainv, Phasev] = deal(nan(size(freq)));
for ii=1:numel(freq),
    fr = freq(ii);
    fr = clip(fr, min(ffreq(:)), max(ffreq(:))); % avoid out-of-range errors 
    jfr = f2j(fr);
    spl = SPL(ii);
    gg = fGG(:, jfr); % gain for this tone, using all candidate zwuis SPLs 
    pp = fPP(:, jfr); % phase for this tone, using all candidate zwuis SPLs 
    Pxout_db = 30+spl+gg-P2dB(fr); % displ power
    [Resx(ii), ibest] = min(abs(Pxout_db-fPx_dB(:,jfr)));
    SPLeqx(ii) = fSPL(ibest);
    Gainx(ii) = gg(ibest);
    Phasex(ii) = pp(ibest);
    %==
    Pvout_db = spl+gg; % vel. power
    [Resv(ii), ibest] = min(abs(Pvout_db-fPv_dB(:,jfr)));
    SPLeqv(ii) = fSPL(ibest);
    Gainv(ii) = gg(ibest);
    Phasev(ii) = pp(ibest);
end
Phasex = cunwrap(Phasex);
Phasev = cunwrap(Phasev);
Phasex = Phasex-repmat(round(Phasex(1,:)),numel(Phasex(:,1)),1);
Phasev = Phasev-repmat(round(Phasev(1,:)),numel(Phasev(:,1)),1);

S = CollectInStruct('-tonestim', freq, SPL, Rcrit, ...
    '-displacement', SPLeqx, Gainx, Phasex, Resx, ...
    '-velocity', SPLeqv, Gainv, Phasev, Resv);

Gain_zw = cat(1,A.Gain) + pmask(cat(1,A.Alpha)<=Rcrit);
Phase_zw = delayPhase(cat(1,A.Phase) + pmask(cat(1,A.Alpha)<=Rcrit),cat(1,A.Fprim)./1e3,tau);
SPLzw = cat(1,A.Gain) + pmask(cat(1,A.Alpha)<=Rcrit);
Z = CollectInStruct(Fprim, Gain_zw, Phase_zw, SPLzw);

