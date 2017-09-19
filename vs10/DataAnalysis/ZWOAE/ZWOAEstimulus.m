function [Z1, Z2, StimParam] = ZWOAEstimulus(igerbil, idataset);
% ZWOAEstimulus - stimulus waveforms of ZWOAE dataset
%    [Z1, Z2] = ZWOAEstimulus(igerbil, idataset) returns the complex
%    analytic stimulus waveform of the two primary groups corresponding to
%    F1 and F2 that was used in ZWOAEdataset (igerbil,idataset). 
%    Only a single period of the stimulus is computed.
%    By convention, 0 dB SPL corresponds to RMS = 1. The real-valued
%    waveforms are obtained by taking the real parts of Z1 and Z2.
%
%    [Z1, Z2, StimParam] = ZWOAEstimulus(..) also returns the stimulus
%    parameters in a struct StimParam containing the following fields:
%        dt: sample period [ms] of the waveforms
%      Fsam: sample rate [kHz] 
%        F1: frequencies [kHz] of the F1 primary tone(s)
%        F2: frequencies [kHz] of the F2 primary tone(s)
%       ph1: starting phases [cycles] of the F1 primary tone(s)
%       ph2: starting phases [cycles] of the F2 primary tone(s)
%        L1: intensity [dB SPL] per tone in F1.
%        L2: intensity [dB SPL] per tone in F2.
%        A1: linear amplitude(s) of primary tone(s) in F1.
%        A2: linear amplitude(s) of primary tone(s) in F1.
%      Nsam: # samples in the stimulus waveforms
%        df: spectral spacing [kHz] of the fft of Z1 and Z2.
%        
%    components. 
%
%    Note: phases are random and do not match those of the original
%    stimulus (don't know how to retrieve them yet).
%
%    See also ZWOAEstruct, getZWOAEdata, REAL.

% extract stimulus params
D = getZWOAEdata(igerbil,idataset, '-nosig');
SP = D.stimpars;
i1 = logical(SP.N1);
i2 = logical(SP.N2);
F1 = SP.F1(i1)/1e3; % Hz -> kHz
F2 = SP.F2(i2)/1e3; % Hz -> kHz
L1 = SP.L1;
L2 = SP.L2;
Fsam = D.fs;  % in kHz
% derived params
A1 = sqrt(2)*dB2A(L1); % srqrt(2) bridges difference pure-tone amplitude vs RMS
A2 = sqrt(2)*dB2A(L2);
A1 = SameSize(A1, F1);
A2 = SameSize(A2, F2);
N1 = length(F1);
N2 = length(F2);
Nsam = D.periodicity;
% random phases
ph1 = rand(size(F1));
ph2 = rand(size(F2));

% generate stim
dt = 1/Fsam; % sample period in ms
t = dt*(0:Nsam-1); % time axis in ms
Z1 = 0;
for ii=1:N1,
    fr = F1(ii); 
    ph = ph1(ii);
    A = A1(ii);
    Z1 = Z1 + A*exp(2*pi*i*(fr*t+ph));
end
Z2 = 0;
for ii=1:N2,
    fr = F2(ii); 
    ph = ph2(ii);
    A = A2(ii);
    Z2 = Z2 + A*exp(2*pi*i*(fr*t+ph));
end

df = Fsam/Nsam;

StimParam = CollectInStruct(dt, Fsam, F1, F2, ph1, ph2, L1, L2, A1, A2, Nsam, df);


