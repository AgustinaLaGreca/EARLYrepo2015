function Ph= Phase(T);
% transfer/phase - phase of transfer object
%    phase(T) returns an array containing the unwrapped phases (in cycles)
%    of the transfer function in transfer object T. 
%
%    See Transfer, Transfer/measure, Transfer/frequency.

if ~isfilled(T),
    error('Transfer object T is not filled.');
end

Ph = cunwrap(cangle(T.Ztrf));



