function Ab = actabit(figh, Name, varargin)
% actabit - construct small actor
%    Actabit(figh, 'Joe', 'wrapup', @foo) returns an Action object 
%    named Joe stored in graphics object with handle figh, which calls 
%    function foo at wrapup time. That's all it does. Actabit can be viewed 
%    as a kind of anonymous class; it dispenses with the need to create an 
%    entire new Action subclass that does a trivial thing, which is a 
%    hassle because you have to tell what *not* to do all the time. 
%  
%    Actabit(figh, 'Joe', 'wrapup', {@foo arg1 ..}) passes arguments to
%    function foo, i.e. it calls foo(arg1, ..). Note that the actibit
%    object itself never enters the story. The callbacks like @foo do not
%    have access to the Actabit object, which really is a dummy.
%   
%    Actabit(figh, 'Joe', 'prepare', {@foo arg1 ..}, 'clear', @fee, ...) 
%    calls foo(arg1,..) at prepare time; fee() at clear time, etc.
%
%    Note that Actabit does not launch a timer at start time (unless you 
%    explicitly tell it to do so :). In that sense it is a rather trivial
%    Action object. Its status is immediately set to 'finished' at start
%    time. When stopped, however, its status is set to 'stopped' anyhow.
%
%    The possible Stages are (quotes omitted): 
%        prepare, start, stop, wrapup, wrapupstopped, clear.
%
%    Here, 'wrapupstopped' overrides any wrapup call when the Actabit has
%    status 'stopped' at the time of the wrapup call. This allows you to
%    differentiate between 'finished' and 'stopped' cases. If no wrapupstop
%    call is defined, Actabit/wrapup will simply call the specified wrapup
%    call (if any), regardless of the status.
%    
%   
%    See also Action, Action/prepare.

% store the desired action in a struct. Field names are stages, field values are actions.

AllStages = {'prepare' 'start', 'stop', 'wrapup' 'clear'};

Ab.Todolist = struct([]); % if varargin is empty, Ab is still well defined
for itrail=1:2:nargin-2,
    stage = lower(varargin{itrail});
    if isempty(strmatch(stage, AllStages, 'exact')),
        error(['Invalid stage ''' stage '''.']);
    end
    Ab.Todolist(1).(stage) = cellify(varargin{itrail+1});
end
Ab = class(Ab, mfilename, action('initialized',nan));
Ab = upload(Ab, figh, Name);





