function prepare(G);
% grabevents/prepare - prepare Grabevents object for its job
%    prepare(G) prepares event grabbing by Grabevent object G.
%    The preparation consists of the following actions:
%        1. check to make sure that G's status is 'initialized'.
%        2. check to make sure that no SeqPlay is currently active;
%        3. perform a calibration of the hardware timing offsets;
%        4. if the source device has gettable settings, get them.
%        5. enable event grabbing by the Seqplay circuit
%        6. reset the event buffers
%        7. set G's status to 'prepared'.
%
%   As with most Action objects, the real action is handled by the Oneshot
%   method, and the timer for spawning calls to Oneshot is created only at
%   start time.
%
%   See also grabevents/start grabevents/oneshot, seqplaygrab.

eval(IamAt);

% 1. check to make sure that G's status is 'initialized'.
[G, Mess] = criticalDownload(G, 'initialized', 'prepare');
error(Mess);

% 2. check to make sure that no SeqPlay is currently active;
SPstat = seqplaystatus;
if SPstat.Active, 
    error('Cannot prepare event grabbing while SeqPlay is active.'); 
end
% 3. perform a calibration of the hardware timing offsets;
G = setdatabuf(G, 'TimingCalib', TimingCalibrate(G)); % determine time lag between DAC and event timer

% 4. if the source device has gettable settings, get them.
G = getsourcedevicesettings(G);
upload(G); % anticipate downloading by status call below

% 5. enable event grabbing by the Seqplay circuit
seqplaygrab on;

% 6. reset the event buffers
feval(G.Hardware.ResetStamping_trigger); % see recHardware

% 7. set E's grabstatus to 'prepared'.
status(G, 'prepared'); % uploading implicit in this call

