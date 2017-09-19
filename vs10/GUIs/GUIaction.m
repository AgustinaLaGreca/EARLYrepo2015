function GUIaction(h, fcn, varargin);
% GUIaction - delegate work to action objects of GUI
%    GUIaction(figh, @foo) finds all action object stored in GUI figh and
%    tells them to foo. The action objects are retrieved using
%    GUIactionList and are visited one by one (in the order they have been
%    uploaded to the GUI).
%
%    GUIaction(figh, @foo, @wrapup, DT) starts a timer that checks every
%    DT ms whether the foo-ing action objects have stopped or finished.
%    Once they all have, the timer tells them to wrapup, one by one.
%    Default DT=80 ms. 
%
%    GUIaction(figh, @foo, @wrapup, @clear, ..., DT) also invokes @clear, ...
%    on all listed Action objects, one by one. Note that all objects are
%    first visited by @wrapup before they are visted by @clear, etc.
%
%    GUIaction(figh, @foo, {@wrapup arg1 arg2 ..}, ...) uses @wrapup with
%    addition arguments, i.e., when the time is there, the call is
%    wrapup(A, arg1, arg2, ..)
%
%    See GUIactionList, action, action/isready, action/isstopped.


[L, N] = GUIactionList(h);
eval(IamAt);


% apply fcn to list members in order
local_apply_to_list(fcn, L, N);

if nargin>2,
    if isnumeric(varargin{end}), % GUIaction(h, fcn, @A, ...@Z, dt_wrp)
        FinalFcns = varargin(1:end-1);
        dt_wrp = varargin{end};
    else,   % GUIaction(h, fcn, @A, ...@Z)
        FinalFcns = varargin;
        dt_wrp = 80; % ms default check-for wrapup interval
    end
    h = timer('Name', 'Jack the Wrapper', ...
        'Period', dt_wrp/1e3, ...
        'ExecutionMode', 'fixedSpacing', ...
        'TasksToExecute', inf, 'StopFcn', '');
    % now that we know h, we can pass it to the timer itself, so it can
    % self-destruct
    set(h, 'TimerFcn', {@local_wrp_callback, FinalFcns, L, N, h});
    start(h);
end

%===============================
function local_wrp_callback(Src, Event, FinalFcns, L, N, h); %#ok<INUSL>
% wrap-up calls caught in a single callback
if local_allSilent(L,N),
    eval(IamAt);
    stop(h); delete(h); % close the door after you (h is handle of timer invoking this very call)
    for ii=1:numel(FinalFcns), % use functions one by one, to visit all listed actions in order
        local_apply_to_list(FinalFcns{ii}, L, N);
    end
end
        
function local_apply_to_list(fcn, L, N);
% apply fcn to members from list one by one
for ii=1:N, % call the individual action
    Ax = L(ii).address;
    if iscell(fcn), % {@fcn arg1 arg2 ..}
        fcn{1}(get(Ax), fcn{2:end});
    else,
        fcn(get(Ax));
    end
end

function as = local_allSilent(L,N);
% checks whether all action by L-members is over
% note: must either *all* be finished, or *all* be stopped!
% The idea is that the stop action must complete visiting all listed 
% actions before we give permission ("all silent - go ahead") to kick in.
% Some actions, like actabit, are finished once they're started. If a mixed
% bag of finished & started is enough to trigger an "all silent", it may
% will interrupt ongoing stopping action. In that case, the stopping action
% must be given time to also stop the 'finished' actions and change their
% status to 'stopped'. The alternative is to introduce something like an 
% intermediate 'stopping' status, but that is quite contrived.
as = 0;
all_finished = 1;
all_stopped = 1;
for ii=1:N, % call the individual action
    Ax = L(ii).address;
    St = status(get(Ax));
    all_finished = all_finished && isequal('finished', St);
    all_stopped = all_stopped && isequal('stopped', St);
    if ~all_finished && ~all_stopped, return; end % it's not getting better
end
% if we get here, all are silent
as = 1;



