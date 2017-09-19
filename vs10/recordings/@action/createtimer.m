function t = createtimer(A)
% createtimer - create timer for action object
%    h = createtimer(A) creates a timer with handle h for action object A.
%    Its TimerFcn will be a tandem of oneshot(A) followed by 
%    finishWhenReady(A). In order for this to work if a given subclass, 
%    both oneshot and isready must be overloaded (the latter is called by 
%    finishWhenReady). Its properties can be accessed via action/thandle 
%    and, implicitly, by action/DT.
%
%    Createtimer is typically invoked by action/start.
%
%    See action, playsound/prepare.

%
%    Note: ugly implementation, but I see no alternative. Creating the
%    timer @ constructor time is no option, because A has to be passed to
%    the callback in its full subclass glory - not just its generic Action part.
%    On the other hand, the timer properly belongs to the Action part of A.
t = timer('Name', class(A), ...
    'TimerFcn', {@oneshotplus A}, 'Period', A.DT(1)/1e3, ...
    'Startdelay', A.DT(2)/1e3, 'ExecutionMode', 'fixedSpacing', ...
    'TasksToExecute', inf, 'StopFcn', '');

function oneshotplus(Src, Event, A) %#ok<INUSL>
% tandem as in help text, Ignore obligatory callback args
download(A);
oneshot(A);
FinishWhenReady(A);










