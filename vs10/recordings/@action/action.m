function A=action(Status, DT);
% action - constructor of action object
%    A = action('Foo',DT) returns an action object having status 
%    Foo and, an inter-action interval DT [ms]. An 1x2 array, DT=[dt t0], 
%    specifies both the timer interval dt and the start delay t0 (see TIMER).
%    The default start delay t0 is a random number between 1 and 100 ms.
%    (The idea is that a spread of true start times will prevent clogging.)
%
%    Action is an abstract class, acquiring its usefulness through
%    subclasses like datagrabber, playaudio, actabit, etc.
%    Any action object must at least have the following methods:
%        prepare: prepare the action. This usually does not mean building
%                 it from scratch, but rather adjust hardware settings,
%                 upload samples, etc.
%        oneshot: perform a single action
%          start: start intermittent action in the background (oneshot per DT)
%           stop: interrupt ongoing action
%         resume: resume previously interrupted action
%        isready: tests whether the action is finished
%    FinishWhenReady: test if ready; if so, stop timer & and update status
%         wrapup: cleanup after finishing. This may include saving data,
%                 clearing variables from memory, reporting results, etc.
%         status: set/get status (one of initialized, prepared, 
%                 started, stopped, finished, wrappedup)
%             DT: get DT value. 
%        thandle: get/set timer handle.
%        
%    Two of these functions: oneshot and isready, must be overloaded.
%    for each concrete subclass (isready is a special case). 
%    The other ones have generic versions, which in many cases will need 
%    overloading anyhow.
%
%    A timer will be later created and attached to A. See
%    action/createtimer. Its TimerFcn will be a tandem of 
%    oneshot(A) followed by finishWhenReady(A). In order for this to work 
%    for a given subclass, both oneshot and isready must be overloaded
%    (the latter is called by finishWhenReady). Its properties can be
%    accessed via action/thandle and, implicitly, by action/DT.
%
%    Action is a subclass of Dynamic, so it has the download and upload
%    methods available for pointer-like programming constructs.
%
%    See also TIMER, collectdata, playaudio, dynamic, GUIaccess.

if nargin<1,
    [Status, DT, AfterYou, thandle] = deal([]);
else,
    if numel(DT)==1, DT = [DT RandomInt(100)]; end % default start delay is random nmber between 0 and 100 ms.
    % create timer and attach it to A (see help text)
    thandle = [];
    AfterYou = [];
end


% NOTE: cannot upload @ construction time, because uploaded A will have
% superclass "Action" - not the concrete subclass it should have
S = CollectInStruct(Status, DT, AfterYou, thandle);
A = class(S,mfilename, dynamic);



