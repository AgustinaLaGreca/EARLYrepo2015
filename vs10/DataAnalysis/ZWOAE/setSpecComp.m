function Spec = setSpecComp(df, Spec, Freq, X)
% setSpecComp - set components of spectrum @ requested frequencies
%
%   Spec = setSpecComp(df, Spec, Freq, X) sets the the values in spectrum Spec 
%   at frequencies Freq(k) to X(k). df si the frequency spacing of Spec.
%   Abs(Freq) may not exceed the sample frequency; negative frequencies
%   are okay - they are taken from the upper half of the spectrum. Note that
%   the use of negative frequencies requires that you pass the entire
%   spectrum as returned by fft.
%
%   See also SpecComp, ZWOAEsubspec.

N = numel(Spec);
Fsam = df*N;
if abs(Freq)>Fsam,
    error('Absolute values of frequencies may not exceed sample rate.');
end

iFreq = 1+round(Freq/df);  % indices of components in spectrum
iFreq = 1+mod(iFreq-1,N);
Spec(iFreq) = X;




