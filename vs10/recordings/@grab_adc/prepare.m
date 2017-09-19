function prepare(G);
% grab_adc/prepare - prepare grab_adc object for its job
%    prepare(G) prepares downloading sampels Grab_adc object G.
%    The preparation consists of the following actions:
%        1. check to make sure that G's status is 'initialized'.
%        2. check to make sure that no SeqPlay is currently active;
%        3. perform a calibration of the hardware timing offsets;
%        4. enable sample downloading by the Seqplay circuit
%        5. reset the ADC buffer
%        6. determine the size of the ADC buffer
%        7. get settings of Source Device if needed.
%        8. set G's status to 'prepared'.
%
%   As with most Action objects, the real action is handled by the Oneshot
%   method, and the timer for spawning calls to Oneshot is created only at
%   start time.
%
%   See also grab_adc/start grab_adc/oneshot, seqplaygrab.

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
G = setdatabuf(G, 'TimingCalib', TimingCalibrate(G)); % determine time lag between DAC and ADC

% 4. enable event grabbing by the Seqplay circuit
seqplaygrab_adc on;

% 5. reset the ADC buffers
feval(G.Hardware.ResetADC_trigger); % this action is defined in recHardware

% 6. determine the size of the ADC buffer
G.BufSize = sys3ParTag(G.Dev, G.BufTag, 'TagSize');

% 7. if the source device has gettable settings, get them.
G = getsourcedevicesettings(G);
upload(G); % anticipate downloading by status call below.

% 8. set E's grabstatus to 'prepared'.
status(G, 'prepared'); % uploading implicit in this call









