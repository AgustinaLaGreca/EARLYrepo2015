function start(R);
% playsound/start - kickstart of playsound action
%    start(P) triggers the D/A by giving a go to the hardware
%
%    See also action/stop, playsound/prepare, playsound/start.

eval(IamAt('indent')); % indent action method below
seqplaygo;

% wait for hardware to start, otherwise finishWhenReady (which is implicitly
% launched by start(R.action) might think we're ready & immediately end the play.
tic;
while toc<0.05, % time out needed for D/A so short that we miss it from the software
    St = seqplaystatus;
    if St.Active, break; end
end % .. can only be set once we're sure we have started ;)

start(R.action); % generic starting: provide & start timer, update status



