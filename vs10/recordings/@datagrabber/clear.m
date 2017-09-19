function varargout = clear(G, varargin);
% datagrabber/clear - clear datagrabber object and its data
%   cleardata(G) offloads G and clears from memory any data collected by 
%   datagrabber G. 
%
%   For each concrete Datagrabber subclass, Clear must be overloaded.
%   It is an error to call the non-overloaded version.
%
%   See also Grabevents/clear, action/clear.

error(['Clear is not yet overloaded for subclass ''' class(G) '''.' char(10) ...
    'It is an error to call the non-overloaded datagrabber/clear.']);






