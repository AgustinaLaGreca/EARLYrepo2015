function varargout = stop(G, varargout);
% datagrabber/stop - stop datagrabber object
%   stop(G) stops data grabbing by G
%
%   For each concrete Datagrabber subclass, Stop must be overloaded.
%   It is an error to call the non-overloaded version.
%
%   See also Grabevents/stop, action/stop.

error(['Stop is not yet overloaded for subclass ''' class(G) '''.' char(10) ...
    'It is an error to call the non-overloaded datagrabber/stop.']);






