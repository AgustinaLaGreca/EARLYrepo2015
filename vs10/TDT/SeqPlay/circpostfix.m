function Postfix = circpostfix(Circuit);
% circpostfix - extract postfix from (full) circuit name
%   S = circpostfix(Circuit) extracts the postfix from circuit name
%   Circuit
%
%   See also SeqPlayInit.

% [path, cname, extension, version] = fileparts(Circuit);
[path, cname, extension] = fileparts(Circuit);
sepi = find(cname=='-',1,'first');
Postfix = cname(sepi:end);






