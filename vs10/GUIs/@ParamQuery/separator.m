function Y = separator(Q, Str);
% ParamQuery/separator - get or set separator property of Paramquery object
%   separator(Q) returns the separator string of Q, which is '' when none
%   was set. Q = separator(Q, S) sets the separator string to S, which must
%   be a valid Matlab identifier.
%
%   The separator property of a Paramquery Q is set true by GUIpiece/add when
%   Q is the first query added to the GUIpiece. GUIval uses the separator
%   property to partition the struct that it returns.
%
%   See ParamQuery, GUIval, ISVARNAME, GUIpiece/add.

if ~isequal(1, numel(Q)),
    error('Separator needs scalar paramquery as input.');
end
if nargin<2, % get
    Y = Q.Separator;
    if isempty(Y), Y = ''; end % stringify
else, % set
    if nargout<1, error('Need an output arg when setting separator value'); end
    if ~isvarname(Str),
        error('Second input arg must be valid Matlab identifier.');
    end
    Q.Separator = Str;
    Y = Q;
end



