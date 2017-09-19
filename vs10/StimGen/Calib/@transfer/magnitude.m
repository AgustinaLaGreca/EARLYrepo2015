function M= magnitude(T);
% transfer/level - magnitude of complex transfer function of Transfer object
%    magnitude(T) returns an array containing the magnitude in dB of the
%    complex transfer function of Transfer object T.
%
%    See Transfer, Transfer/measure, Transfer/Frequency, Transfer/Phase.


if ~isfilled(T),
    error('Transfer object T is not filled.');
end

M = A2dB(abs(T.Ztrf));




