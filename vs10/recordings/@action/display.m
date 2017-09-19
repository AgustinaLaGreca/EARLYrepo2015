function display(T)
% Action/display - DISPLAY for Action objects.
%
%   See Dataset.

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    disp(T);
else
    disp(' ');
    disp([inputname(1) ' =']);
    disp(' ');
    disp(T);
end