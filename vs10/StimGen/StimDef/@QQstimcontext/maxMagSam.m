function MX=maxMagSam(C);
% stimcontext/maxMagSam - maximum magnitude of samples for D/A conversion
%     maxMagSam(S) returns the maximum absolute value (magnitude) of the
%     samples to be sent to the DA converter.
%
%     See stimcontext.


switch C.Hardware.DAC,
    case 'RX6', % sample amplitude is Voltage; max abs is 10 V
        MX = 10;
    otherwise,
        error(['Unknown DAC ''' C.Hardware.DAC '''.']);
end
