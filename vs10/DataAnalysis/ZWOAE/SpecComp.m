function [S, iFreq, isNeg, iposFreq] = SpecComp(df, Spec, Freq)
% SpecComp - return components of spectrum @ requested frequencies
%
%   S = SpecComp(df, Spec, Freq) returns the values in spectrum Spec @
%   frequencies Freq. df gives frequency spacing within Spec
%   Abs(Freq) may not exceed the sample frequency; negative frequencies
%   are okay - they are taken from the upper half of the spectrum. Note that
%   the use of negative frequencies requires that you pass the entire
%   spectrum as returned by fft.
%
%   [S, iFreq, isNeg, IposFreq] = SpecComp(df, Spec, Freq) also returns the 
%   indices of the frequencies in Spec, a logical array isNeg indicating which 
%   components of S correspond to negative frequencies, and the indices iposFreq
%   of the positive frequencies corresponding to all the components in S.
%
%   See also ZWOAspec, ZWOAEsubspec, setSpecComp.

N = numel(Spec);
Fsam = df*N;
if abs(Freq)>Fsam,
    error('Absolute values of frequencies may not exceed sample rate.');
end

iFreq = 1+round(Freq/df);  % indices of components in spectrum
isNeg = iFreq<0;
iFreq = 1+mod(iFreq-1,N);
S = Spec(iFreq);

iposFreq = min(iFreq,N+2-iFreq); % 2 due to matlab base-1 array indices




