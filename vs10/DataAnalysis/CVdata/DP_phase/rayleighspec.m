function [r, pval] = rayleighspec(X, dt, Freq, rampdur, Nchunk)
%rayleighspec -return vector strength & conf. level for spectral comps. based on Rayleigh
%
% == SYNTAX ==
%   [r, pval] = rayleighspec(X, dt, Freq, rampdur, Nchunk)
% == INPUTS ==
%   X:        time signal [-]
%   dt:       time steps in signal X [in ms]
%   Freq:     frequency components for which to calculate alpha [in Hz]
%   rampdur:  duration of applied ramps on window [in ms; see simplegate.m)
%   Nchunk:   number of blocks to divide X in when calculating alpha (see
%             RaySign.m)
% == OUTPUTS ==
%   r:       array of complex mean vectors (abs(r) gives vector strength)   
%   pval:    array of min. significance levels based on Rayleigh's z. 
%
% NOTE: X can be an 1D array only.
%
% See also simplegate, RaySign

if ~isvector(X),
    error('Input argument X miust be a 1D array.');
end
flipflag=~(size(Freq,2)>1);             % is Freq a row (0) or not (1) ? 
X = X(:).';                             % force X into a row
Nsam = numel(X);                        % total # samples
NsamChunk = floor(Nsam/Nchunk);         % # samples in one chunk
NsamRamp = round(rampdur/dt);           % # samples for ramps in simplegate.m
chunk_offset = (0:Nchunk-1)*NsamChunk;  % offsets of chunks

df = 1e3/(Nsam*dt);         % freq spacing in Hz
iFreq = 1+round(Freq/df);   % index in fourier spectrum of Freq component(s)

% restrict X to each of the chunks, one by one
alpha = zeros(Nchunk, numel(iFreq));
for ichunk=1:Nchunk,
    irange = chunk_offset(ichunk)+(1:NsamChunk);
    ch = 0*X;
    ch(irange) = simplegate(X(irange), NsamRamp);
    SP = fft(ch);                       % chunk spectrum
    alpha(ichunk,:)=angle(SP(iFreq));  % get phase of Freq components in rad
end

%use angles in alpha to calculate mean vector (r; complex),
r = mean(exp(1i*alpha));
% pval = RaySign(r,Nchunk);   %get the confidence levels

if flipflag, %make outputs same dimension as Freq
    r=r(:);
    pval=pval(:);
end    



