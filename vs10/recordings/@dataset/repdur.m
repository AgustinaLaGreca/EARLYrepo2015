function d=repdur(D, iCond);
% Dataset/repdur - duration of one stimulus presentation
%   repdur(D) returns the durations in ms of all the single stimulus
%   conditions in order of presentation, including any trailing silence 
%   (i.e., the onset-to-onset intervals).
%
%   repdur(D, iCond) returns a single number quantifying the rep dur of
%   stimulus condition iCond.
%
%   Note that any baseline recordings, whether pre- or poststimulus, are 
%   excluded.
%
%   See also Dataset/burstdur.

iCond = arginDefaults('iCond',nan);

SP = D.Stim.Presentation;
d = SP.PresDur;

if ~isnan(iCond), % single condition requested
    if ~isscalar(iCond), error('Input arg iCond must be single number'); end
    d = unique(d(D.Stim.Presentation.iCond==iCond));
    if isempty(d), error('Requested stimulus condition not found in dataset.'); end
    if ~isscalar(d), error('Multiple rep dur values found for requested stimulus condition?!'); end
else, % all presentations except baselines
    % exclude baselines, if present
    if hasbaseline(SP) && SP.PreDur>0,
        d = d(2:end);
    end
    if hasbaseline(SP) && SP.PostDur>0,
        d = d(1:end-1);
    end
end











