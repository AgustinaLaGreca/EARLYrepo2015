function Str = disp(P);
% dataviewparam/disp - disp for dataviewparam objects.
%   disp(P) displays P. 
%
%   Str = disp(P) returns the string that displays P in Str.

if isvoid(P),
    s = 'void dataviewparam object';
elseif isequal([], P.Param),
    s = ['unspecified' char(P.Dataviewer) ' parameters.']
else,
    sp = disp(P.Param);
    s = strvcat(['' char(P.Dataviewer) ' parameters:'], sp);
end

if nargout>1, 
    Str = s; 
else, 
    disp(s);
end





