function N = nsam(T)
% tsig/nsam - number of samples of tsig object.
%
%   nsamT(T) or T.nsam returns the number of sample in tsig object T.
%   If T has N channels, a 1xN vector is returned containing the respective
%   sample counts. By convention , nsam returns NaN for void tsig objects.
%
%   See tsig.

if isvoid(T), N = nan;
else,
    N = cellfun(@length, T.Waveform);
end




