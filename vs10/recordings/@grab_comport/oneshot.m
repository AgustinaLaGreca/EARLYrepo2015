function oneshot(G, varargin);
% grab_comport/oneshot - elementary quantum of action by grab_comport object
%    oneshot(G) performs one call over the serial port and stored the
%    result.
%
%   See also grab_comport/start grab_comport/wrapup, seqplaygrab, action/status.

G = download(G);

if ~isequal('prepared', status(G)) && ~isequal('started', status(G)) ...
        && ~isequal('sloppy', flag),
    error('Status of G must be ''prepared'' or ''started''.');
end

SamBuf = getdatabuf(G.datagrabber, 'Samples'); % hoard to dump the samples
TimeBuf = getdatabuf(G.datagrabber, 'Times'); % hoard to dump the times
% get the values
samval = feval(G.QueryCall{:});
timeval = 1e3*feval(G.TimeCall{:})/G.Fsam; % time in ms
% concatenate it to previous values
cat(SamBuf, samval); 
cat(TimeBuf, timeval); 


