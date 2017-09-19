function display(T)
% comport_data/display - DISPLAY for comport_data objects.
%
%   See comport_data.

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    disp(T);
else
    disp(' ');
    disp([inputname(1) ' =']);
    disp(' ');
    disp(T);
end



