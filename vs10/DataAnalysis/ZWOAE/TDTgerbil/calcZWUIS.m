function [F1,F2,PerSamp,PerDur] = calcZWUIS(N,F1,F2,df,PerDur,fs,n,k)
%calcZWUIS -calculates ZWUIS frequencies for stimulus design
%
% [F1, F2] = calcZWUIS(N, F1, F2, df, PerDur, fs, n)
%
% Input variables:
%
% N      : number of ZWUIS components in F1
% F1     : mean frequency for ZWUIS-group
% F2     : approx. frequency of single primary (Note: F2 can be smaller F1)
% df     : mean spacing between primaries
% PerDur : approx. length of stimulus period [ms]
% fs     : sample frequency [kHz]
% n      : periodicity in the nfoldplusk-function. Defaults to 3
% k      : the remainders modulo n of F1 and F2, respectively: k=[k1,k2]
%
% Output variables:
%
% F1     : ZWUIS frequencies as to be used
% F2     : Single primary frequency as to be used
% PerSamp: # samples for periodicity of stimulus
% PerDur : periodicity of stimulus in ms

if nargin < 7 | isempty(n),
    n = 3;
end
if nargin < 8 | isempty(k),
    k = [1 0];
end

PerSamp = round(PerDur*fs);
PerDur = PerSamp/fs;
Fbase = 1e3/PerDur; %in Hz

df = local_nfoldplusk(Fbase, df, n, 0);
ddf = n*Fbase;
Nhalf = round(N/2);

F1_spacing = df+fliplr(-Nhalf+1:1:(N-Nhalf-1))*ddf;
%F1_spacing = df+(-Nhalf+1:1:(N-Nhalf-1))*ddf;


F1_startatzero = cumsum([0 F1_spacing]);
F1_mean = mean(F1_startatzero);
F1_shift = F1-F1_mean;
F1_shift = local_nfoldplusk(Fbase, F1_shift, n, k(1));
F1 = F1_startatzero + F1_shift;
F2 = local_nfoldplusk(Fbase, F2, n, k(2));

error(ZWOAEtest(F1, F2));

%----------------------------------------------------
function F = local_nfoldplusk(Fbase, Ftarget, n, k)
% rounds Ftarget to nearest integer multiple of Fbase modulo k*Fbase

tmp = round(((Ftarget/Fbase)-k)/n);
F = (tmp*n+k)*Fbase;










