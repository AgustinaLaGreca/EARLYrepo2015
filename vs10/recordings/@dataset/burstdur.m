function d=burstdur(D, iCond)
% Dataset/burstdur - burst duration duration of one stimulus condition
%   burstdur(D, iCond) returns the burst duration in ms of stimulus
%   condition iCond. Note that this may be a pair of numbers [left right].
%
%   See also Dataset/genericStimParams.

if ~isscalar(D), error('Input arg D must be single dataset.'); end
if ~isscalar(iCond), error('Input arg iCond must be single number.'); end

d = GenericStimparams(D,'BurstDur',iCond);
DAC = GenericStimparams(D, 'DAC');
switch DAC(1)
    case 'L', d = d(1);
    case 'R', d = d(2);
end










