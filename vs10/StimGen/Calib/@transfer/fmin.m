function Freq= fmin(T);
% transfer/fmin - min frequency of transfer object
%    fmin(T) returns the min frequency (Hz) used during its measurement.
%
%    See Transfer, Transfer/fmin, Transfer/frequency.


if ~isfilled(T),
    error('Transfer object T is not filled.');
end

Freq = min(T.Freq);




