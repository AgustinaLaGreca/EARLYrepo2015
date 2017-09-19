function [A, dBr] = Resp_amp(T, varargin);
% transfer/Resp_amp - amplitude of response during transfer measurement
%    L = Resp_amp(T) returns an array containing the response levels in 
%    dB re the dBref_resp of T. The elements of L correspond to the elements 
%    of Frequency(T).
%
%    [L, dBr] = Resp_amp(T) also returns the dBref_resp of T.
%
%    See Transfer, Transfer/Stim_amp, Transfer/Frequency.


if ~isfilled(T),
    error('Transfer object T is not filled.');
end

% first reconstruct the *stimulus* amplitude as used by Transfer/measure
Ainfo = T.CalibParam.Amp;
if numel(Ainfo.Amp_stim)==1, % constant output
    InstAmp = Ainfo.Amp_stim;
else, % output amplitude varies with frequency; use interpolation to get instantaneous values
    InstAmp = interp1(Ainfo.Freq, Ainfo.Amp_stim, T.Freq);
end

% compute the response amplitude from the stimulus amp & transfer
A = A2dB(InstAmp.*abs(T.Ztrf));
dBr = T.dBref_resp;




