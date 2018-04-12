function S = wrapup(S);
% savedataprogress/wrapup - one final update of the online data analysis
%    wrapup(S) updates the online analysis once after data collection is
%    complete. This is done by calling oneshot(S, 'force'). After that, a
%    call to action/wrapup makes sure that the timer of S is stopped and
%    its status is set to wrappedup.
%
%    See also savedataprogress/oneshot, action/wrapup.

eval(IamAt('indent')); % indent call to action/wrapup below
oneshot(S, 'force');

% wrapup(S.Action); % stop timer & set S.Status 
wrapup(S.action); % stop timer & set S.Status % .Action -> .action (Jan,
% April 2018)


