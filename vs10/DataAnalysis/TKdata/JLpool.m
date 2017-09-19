function varargout = JLpool(flag);
%  JLpool - database of JL's ABF data
%    S = JLpool returns struct array S containing info that uniquely links
%    Jeannette's ABF files with their SGSR counterparts.
%  
%    For update, list, .. options, see TKpool.
%
%    See also TKpool.

if nargin<1, flag='return'; end

[varargout{1:nargout}] = TKpool(flag, 'JL');

