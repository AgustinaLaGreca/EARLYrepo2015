function T = fliplr(T)
% tsig/fliplr - swap channels of tsig object
%
%   fliplr(T) reverses the order of the channels of tsig object T.
%   This includes their onset.
%
%   See tsig, tsig/flipud, "methodshelp tsig".

T.Waveform = T.Waveform(end:-1:1);
T.t0 = T.t0(end:-1:1);


