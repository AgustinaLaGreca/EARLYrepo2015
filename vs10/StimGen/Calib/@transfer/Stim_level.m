function [L, dBs] = Stim_level(T, varargin);
% transfer/Stim_level - level of stimulus during transfer measurement
%    L = Stim_level(T) returns an array containing the RMS stimulus levels 
%    in dB re the dBref_stim of T. The elements of L correspond to the 
%    elements of Frequency(T).
%
%    [L, dBs] = Stim_level(T) also returns the dBref_stim of T.
%
%    See Transfer, Transfer/Resp_level, Transfer/DA_amp.

if ~isfilled(T),
    error('Transfer object T is not filled.');
end

% first reconstruct the *stimulus* amplitude as used by Transfer/measure
Ainfo = T.CalibParam.Amp;
if numel(Ainfo.Amp_stim)==1, % constant output
    A = Ainfo.Amp_stim + 0*T.Freq; % add 0*freq to impose right size
else, % output amplitude varies with frequency; use interpolation to get instantaneous values
    A = interp1(Ainfo.Freq, Ainfo.Amp_stim, T.Freq);
end

L = A2dB(sqrt(0.5)*abs(A)); % sqrt(1/2) Peak amp -> RMS
dBs = T.dBref_resp;




