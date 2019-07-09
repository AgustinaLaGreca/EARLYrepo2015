function S = zwuisspectrum(P, Freq, Fscale, ADCdelay, Spikes, DSPdelay);
%   zwuisspectrum - spectrum of pre-processed zwuis response
%     S = zwuisspectrum(P, Freq, ADCdelay, Scale), with P the output of
%     zwuisana (called in "prepare" mode) computes the magnitudes and
%     phases of the requested frequency components (Freq in Hz) and  their
%     Rayleigh significance according to anarayleigh in output struct S.
%     Fscale is the linear scaling factor applied to the spectral
%     components, default Fscale = 1e3*1e-2 to convert DAC Voltage to mV
%     while accounting for a 100X gain of the neural amplifier. ADCdelay is
%     the stimulus-to-recording delay [ms] as returned by e.g. 
%     dataset/timelag. Default ADCdelay==0.
%
%     S = zwuisspec(P, Freq, ADCdelay, Scale, SPT, EventDelay), with SPT an
%     array of spike times (in ms) adds the same analysis to the spike
%     train SPT [ms].  Eventdelay is the stimulus-to-event-time delay. If
%     SPT originates from offline spike detection using APsort, this is the
%     same as ADCdelay. This is the default value for EventDelay. For
%     online spikes acquired over the digital input (e.g. using an Accutrig
%     or BAK device) the delay may be obtained from dataset/timelag,
%     specifying the digital channel.
%
%     See also zwuisana, dataset/apple, dataset/timelag, dataset/anadata.
NchunkRayleigh = 10;

[Fscale, ADCdelay, Spikes, DSPdelay] ...
    = ArginDefaults('Fscale, ADCdelay, Spikes, DSPdelay', 1e3*1e-2, 0, 'none', []);

Cdelay = 0;  

[Magn, Phase, Alpha, PM] = local_spec(P.R, P.dt, Freq, P.f2i(Freq), ...
    P.StartPhase, P.RampDur, NchunkRayleigh, Fscale);
% estimate "threshold magnitude" giving significant data
mm = sort(denan(Magn+PM));
if isempty(mm), ThrUpper_dB = deal(nan);
else,
    ThrUpper_dB = mm(1); % min magn that is signif
end
ThrLower_dB = replaceEmpty(max(Magn(isnan(PM))), nan); % max magn that is not signif
Thr_dB = (ThrUpper_dB+ThrLower_dB)/2;

if isequal('none', Spikes),
    [VS, vsPhase, vsAlpha, vsPM] = deal([]);
else, % make TTL-like spike waveform 
    NsamWsp = P.NcycleZwuis*size(P.R,1);
    Wsp = zeros(NsamWsp,1);
    Time = Xaxis(Wsp, P.dt) - P.Predur; % ms time from start of steady state part of zwuis stim
    isp = 1+round(Spikes/P.dt); % indices of spikes in Wsp
    isp = isp(betwixt_incl(isp,1,NsamWsp));
    Wsp(isp) = 1;
    Wsp = sum(reshape(Wsp, [], P.NcycleZwuis),2);
    [VS, vsPhase, vsAlpha, vsPM] = local_spec(P.R, P.dt, Freq, P.f2i(Freq), ...
        P.StartPhase, P.RampDur, NchunkRayleigh, 0); % last input arg indicates it's spike data
end

S = collectInStruct('-apple', ...
    Freq, Fscale, NchunkRayleigh, '-ADC', Thr_dB, ThrLower_dB, ThrUpper_dB, ...
    Cdelay, Magn, Phase, Phase, Alpha, PM, '-spikes', VS, vsPhase, vsAlpha, vsPM);
S = structJoin(P, S);

%===============
function [Magn, Phase, Alpha, PM] = local_spec(R, dt, Freq, icmp, StartPhase, RampDur, Nchunk, ScaleFac);
Zsp = rfft(R).';
Phase = cangle(Zsp(icmp)) - StartPhase; % phase (cycle) re stim 
if isequal(0, ScaleFac), % spike times
    Nsp = sum(R);
    Magn = abs(Zsp(icmp))/Nsp;
    Alpha = rayleigh_pval(Magn, Nsp);
else, % analog data
    Magn = a2db(abs(ScaleFac*Zsp(icmp)));
    Alpha = anarayleigh(R, dt, Freq, RampDur, Nchunk); % rayleigh significance of phase locking to stimulus
end
PM = pmask(Alpha<=0.001);








