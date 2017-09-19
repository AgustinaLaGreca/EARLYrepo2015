function display(P);
% dataviewparam/display - display for dataviewparam objects.
%   disp(P) displays P.
%
%   Str = disp(P) returns the string that displays P in Str.

if isequal(get(0,'FormatSpacing'),'compact'),
    sep = '';
else,
    sep = ' ';
end

disp([inputname(1) ' =']);
disp(sep);
if numel(P)==1, 
    disp(P); 
else,
    disp([inputname(1) ' = ' sizeString(size(P)) ' ' class(P)]);
end
disp(sep);




