function oneshot(S);
% stopguiaction/oneshot - conditional stop of all action of GUI
%    oneshot(S) checks its stop flag. If the flag is True, S first stops
%    its own timer (preventing multiple stop actions) and then stops all 
%    the action associated with its GUI by calling GUIaction - including 
%    its own.
%
%    See also stopguiaction/setflag, GUIaction.

S = download(S);
if S.flag,
    eval(IamAt);
    stop(thandle(S)); % prevent next shot by timer
    GUIaction(S.figh, @stop);
end



