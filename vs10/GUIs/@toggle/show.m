function show(T); %
% toggle/show - show toggle, i.e., update rendering of T
%    toggle(T) updates the rendering of T.
%
%    See also, betoggle, toggle/click.

QBLACK = 0.25*[1 1 1]; 
if T.EnableState, 
    COLOR = T.defaultColor;
else,
    COLOR = QBLACK;
end

% note that setting the userdata to T is essential for next visit
set(T.h, 'String', T.Str0, 'BackgroundColor', COLOR, 'userdata', T);







