function [T, okay]=impose(T,Str,force); 
% toggle/impose - impose value on toggle object
%    T=impose(T,Str) sets the current string of T to value Str, provided it
%    is a valid value (see toggle).
%    
%    [T, okay]=impose(T,Str) also returns a logical value okay, which is 1 if
%    Str was successfully imposed on T, 0 otherwise.
%
%    Note that impose does not render the change - it only changes T.
%    Use toggle/show to render the change.
%
%    See also toggle, toggle/show, toggle/enable.

if nargin<3, force = ''; end;
if nargout<1, error('Too few output args.'); end

if ~ischar(Str),
    error('Str argument must be character string.');
end

if isempty(force)
    okay = ismember(Str, T.StrArray);
else
    okay = 1;
end

if okay,
    T.Str0 = Str;
end





