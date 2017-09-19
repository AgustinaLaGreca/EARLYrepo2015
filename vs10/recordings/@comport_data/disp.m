function Str = disp(D);
% comport_data/disp - DISP for adc_data objects.
%   
%
%   See also adc_data/samples, adc_data/samperiod.

%S = strvcat('0-0-0-0-0 adc_data object 0-0-0-0-0', disp(struct(D)));

RI = recordinstructions(supplier(D));
S = sprintf('comport_data object holding %s data from %s queried each %0f ms.', ...
    RI.DataType, RI.datafieldname, RI.grabInterval(1));

if nargout<1,
    disp(S);
else,
    Str = S;
end














