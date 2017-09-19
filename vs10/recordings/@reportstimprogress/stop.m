function stop(R, MessMode);
% reportstimprogress/stop - stop ongoing reportimng on stimulus progress
%    stop(R) stops reporting the stimulus progress by R.
% 
%    Two steps:
%       1. display 'D/A stopped' to the GUI
%       2. call generic action/stop
%
%    See also sortConditions, GUImessage, stimpresent/DAstatus.

eval(IamAt('indent')); % indent action method below

GUImessage(R.figh, 'D/A interrupted', 'warning');
stop(R.action); % generic action/stop - includes uploading



