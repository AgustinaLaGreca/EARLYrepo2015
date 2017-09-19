function oneshot(G, maxNevent, flag);
% grabevents/oneshot - elementary quantum of action by Grabevents object
%    oneshot(G) grabs the currently available events from the hardware.
%
%    oneshot(G,Nmax) grabs at most Nmax pending events from the hardware.
%
%    oneshot(G, inf, 'sloppy') ignores the test of G's status, which must
%    normally be either 'prepared' or 'started'.
%
%    Typically, oneshot is called implictly by the timer launched by
%    grabevents/start. There is, however, also the "final call" at wrapup
%    time that collects any pending data.
%
%   See also grabevents/start grabevents/wrapup, seqplaygrab, action/status.

if nargin<2, maxNevent = inf; end; % unlimited number of events
if nargin<3, flag=''; end

G = download(G);

if ~isequal('prepared', status(G)) && ~isequal('started', status(G)) ...
        && ~isequal('sloppy', flag),
    error('Status of G must be ''prepared'' or ''started''.');
end

% how many to grab?
et = getdatabuf(G, 'EventTimes'); % this is a hoard, i.e. pointer!
Ngrabbed = numel(et); % # events already uploaded
Ntimed = sys3getpar('EventCount', G.Dev); % total # events timed by device
NtoGrab = min(maxNevent, Ntimed-Ngrabbed); % # events to be grabbed now

% grab'em
Etick = sys3read('EventTick', NtoGrab, G.Dev, Ngrabbed,'I32'); % sample index
Estamp = sys3read('EventTstamp', NtoGrab, G.Dev, Ngrabbed); % offset from start-of-sample in s (contr to TDT doc om TimeStamp compnt)
% finally, concatenate new times to the ones already grabbed
TC = getdatabuf(G, 'TimingCalib');
cat(et, Etick/TC.Fsam + 1e3*Estamp - 1e-3*TC.Lag_us, 2); % event times in ms; correct hardware lag



