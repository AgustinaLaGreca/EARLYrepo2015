function ActionToggleCB(Src, Event, dum, TrueCallback);
% ActionToggleCB - callback for toggling ActionButtons
%
%   ActionToggleCB(src, Event, dum, TrueCallback) is the default "wrapper 
%   callback" for toggling ActionButtons, i.e., ActionButton objects having a
%   cellstring value for their String property.
%
%   ActionToggleCB does the following things in order:
%      1. store the String of the button at click time in the
%         PreviousString property of the ActionButton
%      2. update the string property of the button control
%      3. update the userdata property of the uicontrol to match the
%         changed value of the button.
%      4. call the "true" callback of the ActionButton, that is, the
%         Callback as specified at contruction time
%
%   The true callback has access to both the string value of the button at
%   the time when clicked (A.PreviousString) and the new string value 
%   (A.CurrentString).
%
%   See also ActionButton, ActionButton/draw.

% Src, Event, dum, TrueCallback

% retrieve ActionButton object
A = getGUIdata(Src);
% toggle the string value
A=localToggleString(A); 
setstring(A); % update display and userdata value

% pass args to true callback
if isfhandle(TrueCallback{1}), % insert Matlab style default args Src and Event
    feval(TrueCallback{1}, Src, Event, TrueCallback{2:end});
else, % only pass what was passed
    feval(TrueCallback{:});
end

% update the attached userdata of A

%===================================================
function A=localToggleString(A);
% store cuurrent string, find next string, return updated action button A
A.PreviousString = A.CurrentString;
Nstr = length(A.String);
imatch = strmatch(A.CurrentString, A.String,'exact');
if isempty(imatch), imatch = Nstr; end;
inew = 1+rem(imatch, Nstr); % rotate
A.CurrentString = A.String{inew};

