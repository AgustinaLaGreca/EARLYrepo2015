function display(T)
% StimPresent/display - DISPLAY for StimPresent objects.
%
%   See StimPresent.

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    disp(T);
else
    disp(' ');
    disp([inputname(1) ' =']);
    disp(' ');
    disp(T);
end