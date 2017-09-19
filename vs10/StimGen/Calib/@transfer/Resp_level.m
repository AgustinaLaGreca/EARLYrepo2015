function [L, dBr] = Resp_level(T);
% transfer/Resp_level - level of response during transfer measurement
%    L = Resp_level(T) returns an array containing the response levels in 
%    dB re the dBref_resp of T. The elements of L correspond to the elements 
%    of Frequency(T).
%
%    [L, dBr] = Resp_level(T) also returns the dBref_resp of T.
%
%    See Transfer, Transfer/Stim_level, Transfer/Frequency.


if ~isfilled(T),
    error('Transfer object T is not filled.');
end

% compute the response amplitude from the stimulus level & transfer
L = Stim_level(T) + A2dB(abs(T.Ztrf));
dBr = T.dBref_resp;




