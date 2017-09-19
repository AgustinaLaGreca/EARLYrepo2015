function TDTgerbil(Flow, Fhigh, dtau, Fcutoff_LP);
% TDTgerbil - stuff a gerbil up Tim's ass

Dev = 'RX6';
if nargin<4, Fcutoff_LP = 1000; end; % Hz cutoff of lowpass filter
Fsam = sys3loadCircuit('tauGerbil', Dev, 25);
sys3setpar(Fcutoff_LP, 'Fcutoff_LP', Dev)


if Fhigh>Fsam/2,
    error('Flow exceeds Nyquist rate.');
end
% compute impulse response of All-pass filter
Nallpass = sys3ParTag('','Chirp_IR', 'TagSize');
APdur = Nallpass/Fsam %ms 
IR = [1 zeros(1,Nallpass-1)];
Hf = fft(IR); % complex transfer function = Fourier spectrum of IR
df = Fsam/Nallpass; % spectral spacing of Hf in ms
ph = cunwrap(angle(Hf)/2/pi); % phase spectrum in Cycles

% generate phases that result in a linear increment of group delay between
% .. Flow and Fhigh
tau = zeros(1,Nallpass-1); % start: tau=0 everywhere
ilo = 1+round(Flow/df); ihi = 1+round(Fhigh/df); % indices of Flow and Fhigh
tauslope=linspace(0,dtau,ihi-ilo+1); % linear slope in interval [Flow, Fhigh]
tau(ilo:ihi) = tauslope;
tau(ihi+1:end) = dtau;
tau = tau -mean(tau) + abs(APdur/2);
%dplot(df,tau);
phase = df*cumsum([0 tau]);
%dplot(df,phase);
Cspec = exp(-2*pi*i*phase);
IR = real(ifft(Cspec));
%dsize(IR);
IR = IR/sum(IR);
sys3write('zero', 'Chirp_IR', Dev);
sys3write(IR, 'Chirp_IR', Dev);
% dplot(1/Fsam,IR);
f1_mean = 3000;
f2_mean = 3600;
df_mean = 60;
Ntone = 4;
PerDur = 1000; % approx stimulus periodicity
[F1, F2, PerSamp, PerDur] ...
    = calcZWUIS(Ntone,f1_mean, f2_mean, df_mean, PerDur, Fsam);
Phase = rand(1,Ntone+1);
local_Stim_Specify(0.1*[0*[1 1 1 1] 1], [F1, F2], Phase, Dev);
sys3run



%=================
function local_Stim_Specify(Amp, Freq, Phase, Dev);
N = length(Amp);
for itone=1:5,
    if itone<=N,
        sys3setpar(Amp(itone), ['ToneAmp_' num2str(itone)], Dev);
        sys3setpar(Freq(itone), ['ToneFreq_' num2str(itone)], Dev);
        sys3setpar(360*Phase(itone), ['TonePhase_' num2str(itone)], Dev);
    else,
        sys3setpar(0, ['ToneAmp_' num2str(itone)]);
    end
end









