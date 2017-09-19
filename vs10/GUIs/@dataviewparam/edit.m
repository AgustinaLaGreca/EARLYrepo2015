function Pout = edit(P, kw, varargin);
% dataviewparam/edit - launch GUI for editing a dataviewparam object.
%   P = edit(P) launches a GUI for editing P. This GUI is defined in the
%   Dataset method dataviewer(P).
%
%   Note: when the user "cancels" the GUI, a void dataviewparam is
%   returned.
%
%   See also dotraster.

% ==callbacks: void P and extra args==
if isvoid(P) && nargin>1,
    feval(['local_' kw], varargin{:}); % e.g., local_OK(figh)
    return;
end

% ==calls from command line or function==
if isvoid(P) && nargin==1,
    error('Void dataviewparam objects cannot be edited.');
end

if nargin<2, kw = ' '; end
% Launch the GUI using the newest info from P.Dataviewer
[Template, ParamGUI] = P.Dataviewer(dataset(),'params');
DVname = char(P.Dataviewer);
figh = newGUI([DVname '_param'], ['parameters for ' DVname], {});
SG=GUIpiece([DVname '_param'],[],[0 0],[10 4]);
SG = add(SG, ParamGUI); % the param queries
SG = add(SG, local_Action(figh)); % messenger; OK & Cancel buttons
draw(figh, SG);
set(figh, 'windowstyle', 'modal'); drawnow;
GUImessage(figh, kw);
GUIfill(figh, P.GUIgrab);
setGUIdata(figh, 'DataviewParam', P);
set(figh, 'KeyPressFcn', @(Src, Ev)edit(dataviewparam(),'Key', figh, Ev));
Pout = dataviewparam(); % default: void dataviewparam (see help)
waitfor(figh, 'userdata');
if isSingleHandle(figh), % GUI is still there; get edited P & close
    Pout = getGUIdata(figh, 'DataviewParam');
    close(figh);
end

%=====================================
function Act = local_Action(figh);
% Check/Play/PlayRec dashboard
Act = GUIpanel('Act', '', 'backgroundcolor', 0.75+[0 0.1 0.05]);
MessBox = messenger('@MessBox', 'The problem is nothing ',4, ... % the '@' in the name indicates that ...
    'fontsize', 9, 'fontweight', 'bold'); % MessBox will be the Main Messenger of he GUI
dum = dataviewparam();
OK = ActionButton('OK', 'OK', 'XXXXXX', 'Accept current values and close GUI.', @(Src,Ev,LR)edit(dum,'OK', figh));
Cancel = ActionButton('Cancel', 'Cancel', 'XXXXXX', 'Ignore values and close GUI.', @(Src,Ev,LR)edit(dum,'Cancel', figh));
Act = add(Act,MessBox);
Act = add(Act, OK, 'below', [10 5]);
Act = add(Act, Cancel, 'nextto', [10 0]);
Act = marginalize(Act,[3 5]);

function local_OK(figh);
P = getGUIdata(figh, 'DataviewParam');
param = GUIval(figh);
if isempty(param), % invalid input
    return;
end
Template = P.Dataviewer(dataset(),'params');
param = union(Template, param);
P.Param = structpart(param, fieldnames(Template));
P.GUIgrab = param.GUIgrab;
setGUIdata(figh, 'DataviewParam', P); 
% if P hasn't changed, the waitfor has not yet been woken up. Make sure it
% gets awake.
if isSingleHandle(figh),
    setGUIdata(figh, 'doemp', 'kloemp'); % this surely terminates waitfor in main fcn block
end

function local_Cancel(figh);
% return void
setGUIdata(figh, 'DataviewParam', dataviewparam()); % this terminates waitfor in main fcn block

function local_Key(figh, Event);
switch Event.Key,
    case 'return',
        local_OK(figh);
    case 'escape',
        local_Cancel(figh);
end



