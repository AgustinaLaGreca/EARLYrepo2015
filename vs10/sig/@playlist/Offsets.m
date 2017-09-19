function I = Offsets(P,k)
% playlist/Offsets - sample offsets of waveforms played by playlist object.
%
%   F=Offsets(P) or F=P.Offsets returns the starting indices of the 
%   waveforms. F is a 1xlength(P.iWave) vector indexing the samples of the
%   played waveform that correspond to the start of a waveform.
%   F equals [] when P no list has been specified for P (see playlist/list)
%
%   Offsets(P,K) only returns the offstes of the Kth waveform(s) to be
%   played.
%
%   See playlist/Waveform, playlist/iWave, playlist/iLongWave.

if isempty(P.iplay), % no list specified - return empty
    I = [];
    return;
end
iw = iLongWave(P);
nsam = Nsam(P, iw);
I = cumsum([1 nsam(1:end-1)]);

if nargin>1,
    I = I(k);
end



