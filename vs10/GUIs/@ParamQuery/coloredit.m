function coloredit(Q, Okay, h);
% ParamQuery/coloredit - color edit control of rendered paramquery
%   coloredit(Q, Okay) sets the background color of the edit control of
%   paramquery Q to a color determined logical Okay. If Okay==true, the
%   color is GP.EditOkayColor; if ~Okay, the color is first set to
%   GP.BlinkColor and after a pause of GP.ErrorBlinkDur seconds, it is set
%   to GP.EditErrorColor. GP is the struct retruned by
%   GUIsettings('ParamQuery'). Q must be rendered. If ~Okay and the edit 
%   control is disabled, it will be enabled.
%
%   If Q is an array, each element is set to the desired color. The
%   blinking is as synchronous as possible.
%
%   coloredit(paramquery, Q, H) ignores the first (void) argument and
%   instead gets the handles to the edit controls directly from handle
%   array H.
%
%   See also GUIsettings, ParamQuery/read, GUIval.

persistent GSQ; if isempty(GSQ), GSQ = GUIsettings('ParamQuery'); end

if nargin<3, % get handles from Q
    h = edithandle(Q);
end
if isempty(h), return; end % nothing to do - don't do it

if Okay,
    set(h, 'backgroundcolor', GSQ.EditOkayColor);
else, % first set all edits to blink color, then wait, then set to error color
    set(h, 'backgroundcolor', GSQ.EditBlinkColor);
    pause(GSQ.EditBlinkDur);
    set(h, 'backgroundcolor', GSQ.EditErrorColor, 'enable', 'on');
    drawnow;
end
