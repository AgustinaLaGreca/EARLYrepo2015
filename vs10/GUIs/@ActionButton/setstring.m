function setstring(A,Str)
% ActionButton/setstring - set string property of rendered actionButton
%    setstring(A,'Foo') sets the action button's CurrentString
%    property to Foo and updates the userdata of the pushbutton 
%    uicontrol to match the updated A. 
%
%    setstring(A) is the same as setstring(A, A.CurrentString).
%
%    For "simple" action buttons, whose String property is a character
%    string, the String property is also set to Foo. For toggling Action
%    butttons, the String property is left alone: only the CurrentString
%    property is set.
%    
%   See ActionButton, ActionButton/setxdata, ActionButton/highlight.

if numel(A)>1,
    error('A must be single Action button.');
end

if isempty(A.uiHandles) || ~isUIcontrol(A.uiHandles.Button),
    error('Action button A is not rendered.');
end

if nargin<2,
    Str = A.CurrentString;
end
A.CurrentString = Str;

if ischar(A.String),
    A.String = Str;
end

% update uicontrol rendering of pushbutton
set(A.uiHandles.Button, 'String', Str, 'Userdata', A);
% also update any associated pulldown menu
hm = getFieldOrDefault(A.uiHandles, 'Menu',[]);
if ~isempty(hm),
    set(hm,'Userdata',A);
end

drawnow;
figure(parentfigh(A.uiHandles.Button));









