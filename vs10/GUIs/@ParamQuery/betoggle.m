function h = betoggle(h,StrArray,Str0);
% betoggle - make pushbutton a toggling button
%   betoggle(h,{'String1' 'String2' ..}) makes the pushbutton uicontrol
%   with handle h a toggling button. By clicking the button, its string
%   property will rotate through the set String1, String2 .. . The initial
%   value of the string is String1.
%
%   betoggle(h,{'String1' 'String2' ..}, 'StringX') selects initial value
%   StringX.
%
%   See also toggle/click, UICONTROL, paramquery/draw.

if nargin<3, Str0=''; end % -> first string in StrArray (see toggle)
% create toggle object T
ht = double(h);
T = toggle(ht,StrArray,Str0);
% attach T to button
set(double(h),'userdata',T, 'callback',{@click_toggle T 'Left'},'ButtonDownFcn', {@click_toggle T 'Right'});
% h.UserData = T;
% % h.CreateFcn = @click_toggle;
% h.Callback = {@click_toggle};
%     %ht 'Left'};
% h.ButtonDownFcn = {@click_toggle};
% % ht 'Right'};
% render it
show(T);

