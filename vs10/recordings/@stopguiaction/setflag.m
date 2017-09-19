function setflag(S);
% stopguiaction/setflag - set stop flag to True.
%    setflag(S) sets the stop lag of S to True, triggering its timer
%    to stop all the action associated with GUI S.figh.
%
%    See also stopguiaction/setflag, GUIaction.

eval(IamAt);
try,
    S = download(S);
    %status(S)
    if ~isequal('started', status(S)), return; end % stopper itse;f must be started before it can do the stopping
    S.flag = 1;
    upload(S);
    status(S, 'stopped'); % doing this fast prevents intervention from wrapup spawned by uniform 'finished' status (see GUIaction)
catch, % if processing of a fast 2nd stop interrupt is slow, downloading may fail because ...
    disp('Nice try. Do you play ping pong?')  %  ...the stopper has been cleared following the 1st stop interrupt.
end




