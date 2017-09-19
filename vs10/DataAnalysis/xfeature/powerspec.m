function powerspec(dt,Y,varargin);
% powerspec - plot power spectrum in dB
%    powerspec(dt,Y, ...) plots the power spectrum of time signal Y, wich
%    is sampled using asample period of dt ms. Powerspec ues a hann window
%    prior to computing the FFT of Y. The ellipses ... denote additional
%    arguments passed to xplot. If Y is a matrix, the power spectra of its 
%    columns are plotted.
%
%    See also plot, PWELCH.

if nargin<3, varargin = {'n'}; end % new plot

if isvector(Y) && size(Y,2)>1,
    Y = Y(:);
end

Nsam = size(Y,1);
df = 1/(Nsam*dt);

% columnwise hann win
Ww = hann(Nsam);
% for icol = 1:size(Y,2),
%     Y(:,icol) = Y(:,icol).*Ww;
% end
SP = A2dB(abs(fft(Y)));
xdplot(df,SP, varargin{:});





