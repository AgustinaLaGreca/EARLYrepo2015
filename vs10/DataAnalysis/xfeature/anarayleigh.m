function [Alpha, S] = anarayleigh(X, dt, Freq, RampDur, Nchunk);
% anarayleigh - rayleigh test for phase locking of analog waveform
%   Alpha = anarayleigh(X, dt, Freq, RampDur, Nchunk)
%   returns confidence level Alpha for the Fourier component at frequency
%   Freq Hz of the analog time signal X sampled at dt ms. To determine
%   Alpha, X is divided in Nchunk nonoverlapping chunks by setting all the
%   samples of X outside the chunk to zero, and the Fourier
%   components at Freq Hz are computed for each of these chunks. Their
%   phase values are then subjcted to a Rayleigh test. Fourier transforms
%   are computed after windowing the time signals using a cos^2 window
%   of duration RampDur ms. For arrays Freq, Alpha is an array of the same
%   size. 
%
%   See also FFT, VectorStrength, SimpleGate.

if ~isvector(X),
    error('Input argument X miust be a 1D array.');
end
X = X(:);
Nsam = numel(X); % total # samples
df = 1e3/(Nsam*dt); % freq spacing in Hz
Freq_idx = 1+round(Freq/df); % index in fourier spectrum of Freq component(s)
NsamChunk = floor(Nsam/Nchunk); % # samples in one chunk
NsamRamp = round(RampDur/dt); % # samples of ramps
chunk_offset = (0:Nchunk-1)*NsamChunk; % offsets of chunks
% restrict X to each of the chunks, one by one
for ichunk=1:Nchunk,
    irange = chunk_offset(ichunk)+(1:NsamChunk);
    ch = 0*X;
    ch(irange) = simplegate(X(irange), NsamRamp);
    SP = fft(ch); % chunk spectrum
    Csp(1:numel(Freq), ichunk) = SP(Freq_idx); % sp(i,j) is complex freq comp i of chunk j.
end
Phase = cangle(Csp);
for ifreq=1:numel(Freq),
    [VS(ifreq), Alpha(ifreq)] = vectorstrength(Phase(ifreq,:),1e3);
end
VS = reshape(VS, size(Freq));
Alpha = reshape(Alpha, size(Freq));
Freq = Freq(:);
S = CollectInStruct(Freq, Phase, VS, Alpha, Csp);





