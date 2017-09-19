function w = waveform(T, k)
% tsig/waveform - waveform values of tsig object
%
%   waveform(T) or T.waveform is a cell array containing the sampled 
%   waveform vectors of T. Note that waveforms can be either numerical
%   arrays or playlist objects (see tsig).
%
%   waveform(T,k) or T.waveform{k}, where k is a single index, returns the
%   waveform of channel k in a column array.
%
%   See tsig, "methodshelp tsig".

if nargin<2, 
    cellOut = 1; % return cell array (see help)
    k = 1:nchan(T); % default: all channels
else,
    cellOut = 0; % return numeric array (see help)
    if ~isscalar(k), error('Index k must be scalar.'); end
    error(numericTest(k, 'posint/nonnan/noninf', 'Index k is '));
end

w = T.Waveform(k);

if ~cellOut, w = w{1}; end;



