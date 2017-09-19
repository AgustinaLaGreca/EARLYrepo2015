function report(M, Mess, CLR);
% Messenger/report - report a message through messenger object
%   report(M, Mess, CLR) reports the message Mess through messenger M using 
%   foregroundcolor CLR. Mess can be a string or cellstring. Default CLR 
%   is M.uicontrolProps.ForegroundColor (see GUIsettings).
%
%   M must be a single, rendered, messenger object.
%
%   See also Messenger, Messenger/reset, GUImessage, GUIsettings.

if numel(M)>1,
    error('Input argument M must be single-element Messenger object.');
elseif isempty(M.uiHandles) || ~ishandle(M.uiHandles.Text),
    error(['''' M.Name ''' Messenger is not rendered.']);
end

if nargin<3, CLR = M.uicontrolProps.ForegroundColor; end 

htext = M.uiHandles.Text; % handle to text uicontrol

if isequal(0,Mess), % reset call (see Messenger/reset)
    set(htext, 'String', M.TestLine, M.uicontrolProps); 
elseif isequal(-1,Mess), % blank call (see Messenger/blank)
    set(htext, 'String', ' ', M.uicontrolProps); 
else, % regular call as described in above help text 
    set(htext, 'string', Mess, 'ForegroundColor', CLR);
end;

%show
drawnow;
figure(parentfigh(M.uiHandles.Text));


