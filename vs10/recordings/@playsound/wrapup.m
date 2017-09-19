function wrapup(P, MessMode);
% playsound/wrapup - wrapup playsound action
%    wrapup(P) wraps up the action spawned by Playsound object P.
%    Wrapup performs the following steps:
%        1. set attenuators to their max value
%        2. delete the timer that was launched by wrapup/start
%        3. change P's status to 'wrappedup'
%    Steps 2&3 are dealt with by action/wrapup, which also checks for the
%    proper status of P when invoked.
%
%    See also playsound/start, playsound/clear, action/wrapup, 
%    action/thandle.

eval(IamAt('indent')); % indent action method below
P = download(P);

% 1. set attenuators to max (shut up)
SetAttenuators(P.Experiment, 'max');

% 3. delete timer  &  4. change P's status to 'wrappedup'
wrapup(P.action); % delegate - including upload








