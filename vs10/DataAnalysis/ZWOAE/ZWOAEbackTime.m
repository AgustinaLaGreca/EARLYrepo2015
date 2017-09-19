function S = ZWOAEbackTime(igerbil, irec, DPtype, tauCM2AC, DphAC2);
% ZWOAEbackTime - estimate backward (CM-to-OAE) travel time of DPs
%    ZWOAEbackTime(igerbil, irec, DPtype) estimates from a AC/CM pair of
%    recordings the CM-to-ear canal travel time [us] of the DPs of given type.
%    Input Args:
%      igerbil: gerbil index (see ZWOAEimport)
%         irec: recording index of either the AC(oustical) or the CM
%               measurement. The other one is found by ZWOAEcompanions.
%       DPtype: type of distortion product: 'near', 'far', etc. For a full
%               listing of DP types, see ZWOAEsubspec.
%    ZWOAEbackTime returns the estimated backward travel time in us. If irec
%    is an array, an array of travel time estimates is returned.
%    The result of the analysis is returned in a struct with fields:
%      igerbil,etc: all input args
%          tauBack: CM to acoustic delay [us]
%              Ph0: CM to acoustic phase offset [cycle]
%           DPfreq: DP frequencies [kHz]
%             Dphi: CM to acoustic phase differences per DPfreq [cycle]
%              Dev: deviations of fitted and measured Dphi [cycle] 
%           RMSdev: RMS of Dev.
%       meanDPfreq: mean DP frequency [kHz]
%                W: weight factors of Dphi evaluated from SNR ratios.
%
%    ZWOAEbackTime(igerbil, irec, DPtype, tauCM2AC) compensates for the
%    last part of the travel, i.e., from the base of the cochlea to the
%    microphone. If tauCM2AC is a single number, tauCM2AC [us] is simply
%    subtracted from the estimated backward time. If tauCM2AC has multiple
%    components, it is interpreted as a polynomial describing the
%    compensatory delay [us] as a function of DP frequency [kHz].
%
%    ZWOAEbackTime(igerbil, irec, DPtype, Freq, DphAC2CM) also compensates
%    for the last part of the backward travel, but this time by processing
%    linear phase data. Freq is an array of frequencies [kHz] and DphAC2CM
%    holds the corresponding AC-to-CM phase data [cycles], i.e. Phcm-Phac.
%
%    [tauBack, Dev, DPfreq]=ZWOAEbackTime(..) also returns Dev=RMS(D), 
%    where D are the differences between the phase and the linear fit used 
%    to estimate the group delay, and Freq, the mean frequency [kHz] of the
%    DPs under consideration.
%
%    See also LinPhaseFit, ZWOAEsubspec, ZWOAEimport.
if nargin<4,
    tauCM2AC = 0;
end
if nargin<5,
    DphAC2 = [];
end
% if nargin>4,
%     freq = tauCM2AC;
%     %DphAC2;
%     error('5-input-arg syntax not yet implemented')
% end

if numel(irec)>1, % multiple recordinfg pairs; use recursion
    for ii=1:numel(irec),
        S(ii) = ZWOAEbackTime(igerbil, irec(ii), DPtype, tauCM2AC, DphAC2);
    end
    return;
end

%-----------single irec from here-----------

[irecAC, irecCM] = ZWOAEcompanions(igerbil,irec);
[DPfreq, MagnAC, PhaseAC, LnoiseAC,Wac] = ZWOAEsubspec(igerbil, irecAC, DPtype);
[DPfreq, MagnCM, PhaseCM, LnoiseCM,Wcm] = ZWOAEsubspec(igerbil, irecCM, DPtype);
W = 1./(1./Wac+1./Wcm); W = 100*W/sum(W); % combine weight factors

if nargin == 5 & ~isempty(DphAC2),
    tmp = DPfreq<min(tauCM2AC);
    PhaseAC(tmp) = [];    PhaseCM(tmp) = [];
    DPfreq(tmp) = []; W(tmp) = [];
    PhaseAC = local_DphAC2CM(tauCM2AC, -DphAC2, DPfreq, PhaseAC);
end

Dphi = PhaseAC-PhaseCM; % phase diff between CM and AC

if (numel(W)-sum(isnan(W))) < 3 | isempty(Dphi),
    %no fit possible, skip and return bogus values
    [Ph0, Tau,Dev, RSMdev] = deal(0,0,0,0);
    Dphi = zeros(size(MagnAC));
else
    [Ph0, Tau, Dphi, Dev] = LinPhaseFit(DPfreq, Dphi, [-1 2],W);
    RSMdev = rms(Dev); % individual fit residues -> RMS
end
meanDPfreq = mean(DPfreq);

% apply compensation for "last part" of back travel
if numel(tauCM2AC)>1 & nargin < 5, % polynomial coeff. Evaluate to get compensatory delay
    tauCM2AC = polyval(tauCM2AC, meanDPfreq);
    tauBack = 1e3*Tau-tauCM2AC;
else
    tauBack = 1e3*Tau;
end

S = CollectInStruct(igerbil, irec, DPtype, tauBack, Ph0, DPfreq, Dphi, Dev, RSMdev, meanDPfreq, W, tauCM2AC, DphAC2);
S.PRG = 0;
if ~all(Dphi==0)
    S.PRG = phaseResidue(Dphi, zeros(size(Dphi)), W);
    S.PRG = 100*(1-rms(S.PRG));
end



function PhaseAC = local_DphAC2CM(Fr, Comp, Fdp, PhaseAC)
%compensate extra for "reverse" delay
%Fr and Comp are arrays giving, per frequency, AC-to-CM phase data
%These data tell how much "sooner" the CM occured re. AC.
%for compensation, the AC phase data must be "negatively delayed" as given
%by these phase data.
%
%Fdp is array at which freqs phase accumulation is needed and PhaseAC
%are the uncompensated AC phases at these freqs
%
% PhaseAC is also returned, but compensated...

[dum, iF] = unique(Fr); %make Fr unique, otherwise interp1 doesn't work
Fr = Fr(iF);
Comp = Comp(iF);
Cdp = interp1(Fr, Comp, Fdp,'linear','extrap');
PhaseAC = PhaseAC-Cdp; %compensate AC phase by middle-ear transmission

