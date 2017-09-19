function [Y, Rseed] = scramblePhase(Y, Rseed);
% scramblePhase - scrambl phase of time signal.
%    Y=scramblePhase(Y) performs FFT, randomizes the phase of the complex
%    spectrum and performs IFFT.
%
%    scramblePhase(Y, Rseed) uses randomseed Rseed using SetRandState.
%
%   [Y, Rseed] = scramblePhase(Y) returns the random seed.
%
%    See SetRandState.

if nargin<2, Rseed = []; end

if isempty(Rseed),
    Rseed = SetRandState;
end
SetRandState(Rseed);
IR = isreal(Y);

Y = fft(Y);
Y = Y.*exp(2*pi*i*rand(size(Y)));
Y = ifft(Y);
if IR, Y = real(Y); end


















