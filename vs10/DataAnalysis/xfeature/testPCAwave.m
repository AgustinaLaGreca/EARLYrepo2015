function [W, TT] = testPCAwave(alpha, dt, Ntime);
% testPCAwave - test waveforms for PCA analysis
%   [W, T] = testPCAwave(alpha, dt, Ntime)

if nargin<2, dt = 0; end
if nargin<3, Ntime=200; end

TT = linspace(-1,1,Ntime)';
Win = hann(Ntime);
[a1, a2, T, dt, Win] = SameSize(abs(alpha(1,:)), 2*abs(alpha(2,:)), TT, dt, Win);
T = 12*(T-dt);
W1 = (T+0.25).*exp(-a1.*(T+0.25).^2);
W2 = (T-0.5).*exp(-a2.*(T-0.45).^2);
W = Win.*(W1 - W2 + 0.1*sqrt(W1.^2+W2.^2));







