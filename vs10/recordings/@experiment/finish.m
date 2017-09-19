function finish(E, Cm);
% experiment/finish - finish Experiment and end its role as current one.
%    finish(E) sets the state of E to 'finished' and saves E. This is only
%    allowed if E is the current experiment. Note that this causes E to no
%    longer be the current experiment. The finishing is reported to E's
%    logfile.
%
%    finish(E, 'Foo') appends the text Foo to the logfile report.
%
%    See also Experiment/iscurrent, Experiment/resume, Experiment/state.

if nargin<2, Cm=''; end
if ~iscurrent(E),
    error('Experiment is not the current one; it cannot be finished.');
end

E = load(E);
E = status(E, 'Finished', datestr(now));
E = state(E, 'Finished');
save(E, ['Finished ' status(E, 'Finished') ' ' Cm]);









