function Str = disp(D);
% eventdata/disp - DISP for adc_data objects.
%   
%
%   See also adc_data/samples, adc_data/samperiod.

%S = strvcat('0-0-0-0-0 adc_data object 0-0-0-0-0', disp(struct(D)));

RI = recordinstructions(D); 
S = sprintf('Grabbeddata object holding %s data from %s.', ...
    RI.DataType, RI.datafieldname);

if nargout<1,
    disp(S);
else,
    Str = S;
end














