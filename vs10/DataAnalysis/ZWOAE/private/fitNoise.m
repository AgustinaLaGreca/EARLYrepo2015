function Pnoise = fitNoise(MG, Nloops, df)
%fitNoise -fit Noisefloor with polynomial & return coefficients
%
% Pnoise = fitNoise(MG, Nloops, df), where;
% MG is spectrum, Nloops = Nsam/periodicity, and df = freq. spacing [kHz]
% Function fits a polynomial to the noise floor & returns the
% coefficients of this fit in Pnoise. 
%
% Nsam          : # of samples in recorded signal
% periodicity   : # of samples over which recording is periodic
%
% Nloops determines which freq-bins in the spectrum are noise. The noise
% floor is fitted with a polynomial (order specified in function) and from
% that all noise data points < fitted function are removed. Remaining data
% are again fitted with polynomial,...
% In total, this is repeated 4 times before returning coefficients. 
%
% Syntax example:
%
% D.Pnf = fitNoise(MG, Nsam/D.periodicity);
%

%NOTE: changing the used function to describe noise  also requires changing ZWOAEgetNoise.m

polyorder = 10; %order of the fitted polynomial

warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
halfLength = length(MG)/2;
NoiseFloor = MG(2:Nloops:halfLength); %noise spectrum; 1st half
    
ndf = df*Nloops;
nf_freq = df+ndf*(0:length(NoiseFloor)-1).';

% Fit the noisefloor on a logarithmic x-scale
% Use linspace and interp1 to get more points @ lower freqs
nf_freq = log(nf_freq);
tmp_freq = linspace(min(nf_freq), max(nf_freq), numel(nf_freq));
NoiseFloor = interp1(nf_freq, NoiseFloor, tmp_freq);
nf_freq = tmp_freq;

ibig = 1:numel(NoiseFloor); %1st time, include all points
for ii=1:4,
    nf_freq = nf_freq(ibig);
    NoiseFloor = NoiseFloor(ibig);
    Pnoise = polyfit(nf_freq, NoiseFloor, polyorder);
    NFfit = polyval(Pnoise, nf_freq);
    ibig = find(NoiseFloor>NFfit);
end
warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');



