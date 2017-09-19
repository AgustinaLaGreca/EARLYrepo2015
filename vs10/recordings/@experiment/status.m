function S = status(E, FN, X, varargin);
% experiment/status - get/set status of experiment
%    status(E) returns a struct describing the current status of experiment
%    E, as read from disk. By definition, Status returns [] for void E.
%
%    status(E, 'Foo') returns only field Foo of the status.
%
%    E=status(E, 'Foo', X) sets the value of Foo to X. Setting the status has
%    the side effect of updating the 'Modified' field of the status and of
%    saving the experiment. Setting the status is only allowed for the
%    current experiment. Also, the 'State' parameter may not be set by
%    Status: use Experiment/state to do this.
%
%    E=status(E, 'Foo1', X1, 'Foo2', X2, ... ) sets multiple values.
%
%    See also Experiment/save, Experiment/modified, eExperiment/iscurrent, 
%             Experiment/State, Experiment/touch.


if ~isvoid(E), E = load(E); end % upload to get most recent version
if nargin==1, % return whole Status field
    S = E.Status;
    return
elseif nargin==2, % return selected field of status
    if ~isfield(E.Status, FN),
        error(['Unknown Status parameter ''' FN '''.']);
    end
    S = E.Status.(FN);
    return
end

% set: first test whether allowed
if ~iscurrent(E),
    error('The status of an experiment may only be saved when it is the current experiment.');
end

Nset = (nargin-1)/2;
Ei = E; % keep copy of the original for comparison below
for iset=1:Nset,
    if ~isfield(E.Status, FN),
        error(['Unknown Status parameter ''' FN '''.']);
    elseif isequal('State', FN),
        error('It is not allowed to use Experiment/Status to change an experiment state. Use Experiment/State instead.');
    end
    E.Status.(FN) = X;
    if iset==Nset, break; end
    [FN, X] = deal(varargin{2*iset-1}, varargin{2*iset});
end
% by definition, the StatusModified property has changed
E.Status.StatusModified = datestr(now);
save(E);
[dum, Change] = structcompare(Ei.Status,E.Status);
if ~isempty(Change),
    addtostatuslog(E, strvcat('-------------------------------------', disp(Change)));
end
S = E;








