function A=stop(A, flag);
% action/stop - stop action of action object
%    stop(A) stops the action of A. status(A) must equal 'started', 
%    'stopped', or 'finished'. In the last case, stop(A) leaves A
%    unaffected.
%
%    stop(A,'strict') does stops A when status(A) is 'started', but throws
%    an error when A's status is either 'stopped' or 'finished'.
%
%    action/stop merely starts the timer with handle thandle(A) and makes
%    sets A's status to 'stopped'. If more spefic work is needed for a 
%    certain subclass, the stop method must be overloaded.
%
%    See also action/start, action/FinishWhenReady, collectdata, playaudio.

eval(IamAt);
if nargin<2, flag=''; end
[A, Mess]=criticalDownload(A, {'started' 'stopped' 'finished'}, 'stop');
error(Mess);

if isequal('strict', flag) && ~isequal('started', A.Status),
    error('Status of action object must be ''started'' in order to be strictly stopped.');
end

if isequal(status(A), 'started'), % avoid double stops and stop-after-finish
    stop(A.thandle);
    A=status(A,'stopped'); % implicit uploading here
end


