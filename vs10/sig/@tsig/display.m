function display(T)
% tsig/display - DISPLAY for tsig objects.
%
%   See tsig.

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    disp(T);
else
    disp(' ');
    disp([inputname(1) ' =']);
    disp(' ');
    disp(T);
end






