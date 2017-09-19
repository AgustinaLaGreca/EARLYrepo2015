function [w, S] = multisine(dt, Amp, Dur, CF, BW, N);
%   multisine(dt, Amp, Dur, CF, BW, N)
%     dt: sample period in us
%    Amp: [mV] abs peak value of stimulus
%    Dur: [ms] duration in 
%     CF: center freq in Hz
%     BW: bandwidth in Hz
%      N: # components

Ncyc = 10; % min # cycles in stim

CycDur = Dur/Ncyc % ms duration of sin cycle
meanDF = BW/N; % mean comp spacing
NsamCyc = round(1e3*CycDur/dt);
Nsam = NsamCyc*Ncyc; % ms total dur
df = 1e3/CycDur; % Hz freq spacing

F0 = CF-BW/2; % Hz freq of lowest cmp
% convert frequencies & freq spacings into multiples of 3*df
M0 = round(F0/df/3) 
MmeanDF = round(meanDF/df/3) % mean spacing
MDF0 = MmeanDF -round((N-1)/2) % spacing between lowest components
MDF = MDF0+ (0:N-2);

% back to Freqs
F0 = df*(3*M0+1) % lowest comp
DF = df*(3*MDF+1) % spacings
Freq = cumsum([F0 DF])
mean(Freq);

% generate stimulus
w = 0;
Ph0 = rand(1,N); % starting phase in cycles
time = 1e-6*dt*(0:Nsam-1).'; % time axis in s .. COLUMN vector!!
max(time)
for ii=1:N,
    freq = Freq(ii) % current freq [Hz]
    ph0 = Ph0(ii) % current startin gphase [cycle]
    w = w + cos(2*pi*(ph0 + freq*time));
end
% impose peak value
w = Amp*w/max(abs(w));

% dplot(df,a2db(abs(fft(w))));
% xlim([0 2*max(Freq)]);

df = 1e3/Dur % Hz freq spacing of spectrum

S = CollectInStruct(dt, Amp, Dur, CF, BW, N, Ph0, Freq, ...
    CycDur, F0, MDF, N, NsamCyc, df, DF, Ncyc, Nsam, NsamCyc);






