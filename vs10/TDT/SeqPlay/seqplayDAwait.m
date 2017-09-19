function Mess = SeqplayDAwait(TimeOut);
% SeqplayDAwait - wait for seqplay DA conversion to finish
%   SeqplayDAwait waits until the squenced-play DA conversion, previously
%   triggered by seqplayGo has finished. It immediately returns when no DA
%   conversion was triggered.
%
%   SeqplayDAwait(T) uses a timeout of T ms, after which it returns even
%   when DA conversion has not finished. The default T is 10000 ms, that is,
%   10 seconds.
%
%   Mess = SeqplayDAwait(..) returns a message, which is empty on a regular
%   return, but returns an appropriate text message when timeout has been
%   reached, or when no seqplay DA conversion was triggered in the first
%   place.
%
%   See also Seqplay, SeqplayGo.

tic; % start timer
if nargin<1,
    TimeOut = 10e3; % ms default timeout
end
Mess = ''; % optimistic default
NoSeqMess = 'No seqplay DA conversion triggered.';

PI = private_seqPlayInfo;
if isempty(PI),
    Mess = NoSeqMess;
elseif ~isequal('listed', PI.Status);,
    Mess = NoSeqMess;
end

if ~isempty(Mess), 
    return;
end

Tend = TimeOut*1e-3; % timeout in seconds
while 1,
    sps = seqplaystatus;
    if ~sps.Active, 
        break; 
    elseif toc>Tend,
        Mess = 'Timed out.'
        break;
    end
end



