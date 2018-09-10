function sig=generate_schroeder(f0,C,T,fs,high_freq_limit)
% f0      fundamental frequency                                   Unit: Hertz

% C       speed of frequency shift                                Unit: none
% T       length of the schröder phase tone                       Unit: seconds
% high_freq_limit  Highest frequency in complex                Hz
% fs      sampling frequency                                        Unit: Hertz default fs=48000 Hz

% Output
% sig       vector of the computed Schroeder Phase

%  Carney parameters:  
%  F0 = [50 100 200 400]
%  C = [-1 : 0.25 : +1.0], although more recently I'm shortening this to  C = [-1:  0.5 : 1.0]
%  Used duration of 04 sec
%  Note that fundamental was included here - in some psychophysical studies,
%  the fundamental is excluded.

N = floor(high_freq_limit/f0);   % # of components in complex

%% Timevector
t=0:1/fs:T-1/fs;

%% Sum the fundamental and its harmonics to build the Schröder phase complex
sig=zeros(1,fs*T);

for n=1:N
    
    phi_shift=C*pi*n*(n+1)/N;
    
    sig=sig+sin(2*pi*t*f0*n + phi_shift);
end

end