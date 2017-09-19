function S = wrapup(S);
% showdataprogress/wrapup - one final update of the online data analysis
%    wrapup(S) updates the online analysis once after data collection is
%    complete. This is done by calling oneshot(S, 'force'). After that, a
%    call to action/wrapup makes sure that the timer of S is stopped and
%    its status is set to wrappedup.
%
%    See also Showdataprogress/oneshot, action/wrapup.

eval(IamAt('indent')); % indent call to action/wrapup below
oneshot(S, 'force');
wrapup(S.Action); % stop timer & set S.Status


