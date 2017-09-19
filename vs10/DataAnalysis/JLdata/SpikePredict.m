function P = SpikePredict(Vthr, Vmean, Vstd);
%  SpikePredict - predict firing probability from membrane potential
%     P = SpikePredict(Vthr, Vmean, Vstd) uses membrane potential waveform
%     Vmean and its trial-to-trial standard deviation Vstd to predict the
%     instantaneous firing probability. Inputs are
%         Vthr: firing threshold: either a constant or an array having the
%               same size as Vmean and Vstd.
%        Vmean: array holding instantaneous membrane potential waveform
%         Vstd: instantaneous standard deviation of 
%     Vmean, Vstd, Vthr have the same units (e.g. mV).
%
%     Output P is an array having the same size as Vmean holding the 
%     instantaneous firing probability per time bin.
%
%     For the computation of P, Gaussian statistics is eassumed.


% Normalize threshold re Vmean & Vst
Vthr = (Vthr-Vmean)./Vstd;
P = erfc(Vthr/sqrt(2))/2;




















