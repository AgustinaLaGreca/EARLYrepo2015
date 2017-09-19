function oneshot(G, flag);
% grab_adc/oneshot - elementary quantum of action by grab_adc object
%    oneshot(G) grabs the currently available samples from the hardware.
%
%    oneshot(G, 'sloppy') ignores the test of G's status, which must
%    normally be either 'prepared' or 'started'.
%
%    Typically, oneshot is called implictly by the timer launched by
%    grab_adc/start. There is, however, also the "final call" at wrapup
%    time that collects any pending data.
%
%   See also grab_adc/start grab_adc/wrapup, seqplaygrab, action/status.

persistent LastTime
if isempty(LastTime), LastTime=now; end

if nargin<2, flag=''; end

G = download(G);

if ~isequal('prepared', status(G)) && ~isequal('started', status(G)) ...
        && ~isequal('sloppy', flag),
    error('Status of G must be ''prepared'' or ''started''.');
end

% how many to grab?
sams = getdatabuf(G, 'Samples'); % remember that this is a hoard (pointer)
Ngrabbed = rem(numel(sams), G.BufSize); % # samples already uploaded modulo buf size
Nrec = sys3getpar(G.NsamTag, G.Dev); % total # samples recorded modulo buf size
NtoGrab = Nrec-Ngrabbed;
if NtoGrab>=0, % simple: from where we were to were it got
    offset1 = Ngrabbed;
    Ngrab1 = NtoGrab;
    offset2 = 0;
    Ngrab2 = 0;
else, % recording wrapped around; follow it by jumping back
    NtoGrab = G.BufSize + NtoGrab; % this is the real # of samples to grab
    offset1 = Ngrabbed;
    Ngrab1 = G.BufSize-offset1;
    offset2 = 0;
    Ngrab2 =  NtoGrab - Ngrab1;
end
% disp('-------grabbing & catting----------')
% [offset1, Ngrab1, offset2, Ngrab2]
% disp(datestr(now-LastTime, 'SS:FFF')); LastTime=now;
% grab'em
[sam1, sam2] = deal([]);
if Ngrab1>0, 
    sam1 = sys3read(G.BufTag, Ngrab1, G.Dev, offset1 ,'F32');
end
if Ngrab2>0, 
    sam2 = sys3read(G.BufTag, Ngrab2, G.Dev, offset2 ,'F32');
end
cat(sams, [sam1, sam2],2);  % concatenate new samples to the ones already grabbed









