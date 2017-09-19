function prepare(G);
% grab_comport/prepare - prepare grab_comport object for its job
%    prepare(G) prepares event grabbing by Grabevent object G.
%    The preparation consists of the following actions:
%        1. check to make sure that G's status is 'initialized'.
%        2. perform the PrepCall of the device
%        3. set G's status to 'prepared'.
%
%   As with most Action objects, the real action is handled by the Oneshot
%   method, and the timer for spawning calls to Oneshot is created only at
%   start time.
%
%   See also grab_comport/start grab_comport/oneshot, seqplaygrab.

eval(IamAt);

% 1. check to make sure that G's status is 'initialized'.
[G, Mess] = criticalDownload(G, 'initialized', 'prepare');
error(Mess);

% 2. if the source device has gettable settings, get them.
feval(G.PrepCall{:});
upload(G); % anticipate downloading by status call below

% 3. set E's grabstatus to 'prepared'.
status(G, 'prepared'); % uploading implicit in this call

