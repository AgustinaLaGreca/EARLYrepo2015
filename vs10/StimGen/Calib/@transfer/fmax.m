function Freq= fmax(T);
% transfer/fmax - max frequency of transfer object
%    fmax(T) returns the max frequency (Hz) used during its measurement.
%
%    See Transfer, Transfer/fmin, Transfer/frequency.


if ~isfilled(T),
    error('Transfer object T is not filled.');
end

Freq = max(T.Freq);




