function ir=isready(A,St);
% action/status - true when action is finished
%    isready(A) returns True when the action spawned by A is finished.
%
%    Generally, the ready-state of A depends on the specifics of A, and
%    action/isready must be overloaded for a concrete subclass. In some
%    cases, however, termination of the A's action may be determined, not 
%    by A's own affairs, but by the termination of other action objects.
%    Such a dependence can be specified by the action/AfterYou method. If no
%    such dependence is defined, action/isready must be overloaded in order
%    for A to be used.
%
%    Note the difference bertween finished and stopped. Finishing means
%    that A is ready with its task, which may be the case because A1, A2,
%    etc are finished. If any of them is stopped, not finished, this does
%    not coujnt as a reason for A to be ready: Stopping is an active
%    process that should visit all active ('started') Action objects.
%    By the same token, when A's own status is 'stopped', isready(A) is
%    false.
%    
%    See also action/afteryou, action/start, action/wrapup.

A = download(A, 'sloppy');
if isempty(A.AfterYou),
    error(['Action object A does not know when to finish. You must either overload '  char(10) ...
        'the isready method for ' upper(class(A)) ' objects, or provide a termination dependence' char(10) ...
        'for this instance of A by using action/AfterYou.']);
end

if isequal('stopped', A.Status), % stopped -> not finished
    ir = 0;
    return;
end

% dependence: A is ready whenever all of the listed objects are finished
ir = 0;
for ii=1:numel(A.AfterYou),
    Stat = status(get(A.AfterYou(ii)));
    if ~isequal('finished', Stat),
        return; % one failure is enough
    end
end
% we only get here if all objects in the AfterYou list are finished
ir = 1;



