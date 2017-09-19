function display(T)
% adc_data/display - DISPLAY for adc_data objects.
%
%   See adc_data.

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    disp(T);
else
    disp(' ');
    disp([inputname(1) ' =']);
    disp(' ');
    disp(T);
end



