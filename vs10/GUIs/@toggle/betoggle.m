function h = betoggle(h,StrArray,Str0);
% % betoggle - make pushbutton a toggling button
% %   betoggle(h,{'String1' 'String2' ..}) makes the pushbutton uicontrol
% %   with handle h a toggling button. By clicking the button, its string
% %   property will rotate through the set String1, String2 .. . The initial
% %   value of the string is String1.
% %
% %   betoggle(h,{'String1' 'String2' ..}, 'StringX') selects initial value
% %   StringX.
% %
% %   See also toggle/click, UICONTROL, paramquery/draw.



if nargin<3, Str0=''; end % -> first string in StrArray (see toggle)
% create toggle object T
T = toggle(h,StrArray,Str0);
T.Str
T.Str0
% attach T to button
set(h,'userdata',T, 'callback',{@click T 'Left'},'ButtonDownFcn', {@click T 'Right'});
% render it
show(T);

function click(h,dum,dumT,MouseButton); %#ok<INUSL>
% toggle/click - click callback of toggle
%    toggle(h,dum,T,'Left') is called by toggle button with handle h when
%    left-clicked. It sets string property of the button to the the next 
%    string from the list (rotating).
%
%    toggle(h,dum,T,'Right') is called by toggle button when right-clicked.
%    It enables/disables the state of the toggle button.

% get current state of toggle, i.e., the object in userdata clicked button
h = double(h);
T = get(h,'userdata');
switch MouseButton,
    case 'Left', % regular click  : rotate to next string
        imatch = strmatch(T.Str0,T.StrArray, 'exact');
        if isempty(imatch), imatch = T.Nstr; end;
        inew = 1+rem(imatch,T.Nstr); % rotate
        T.Str0 = T.StrArray{inew};
        T=enable(T,1);
    case 'Right', % right-click" "enable/disable" button 
        T=enable(T,'swap');
end
% T has changed; re-render it
show(T);


