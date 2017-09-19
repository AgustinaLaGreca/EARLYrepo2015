function stop(R, MessMode);
% playsound/stop - stop playsound action
%    stop(R) interrupts ongoing sound stimulation by R
%    This consists of the following steps:
%        1. stop the timer T.thandle 
%        2. set attenuators to max
%        3. call seqplayhalt
%        4. set the status of R to 'stopped'.
%
%    See also action/stop, playsound/prepare.

eval(IamAt);
R = download(R);
if isequal('finished', status(R)), % just too late; cancel the stop action
    return;
elseif ~isequal('started', status(R)),
    struct(R)
    struct(R.action)
    error('Can only stop object whose status is ''started''.');
end

% 1. stop timer. (Note that doing this first avoids a competitve 'natural' ending of the action spawned by finishWhenRead)
stop(thandle(R));

% 2. attenuators to max (shut up)
SetAttenuators(R.Experiment, 'max');

% 3. stop D/A conversion
seqplayhalt; % really shut up


% 4. change status of R
status(R, 'stopped');






