function Str = disp(D);
% eventdata/disp - DISP for adc_data objects.
%   
%
%   See also adc_data/samples, adc_data/samperiod.

%S = strvcat('0-0-0-0-0 adc_data object 0-0-0-0-0', disp(struct(D)));

RI = recordinstructions(supplier(D));
S = sprintf('Eventdata object holding %s data from %s @ %0.1f kHz.', ...
    RI.DataType, RI.datafieldname, RI.Fsam/1e3);

if nargout<1,
    disp(S);
else,
    Str = S;
end














