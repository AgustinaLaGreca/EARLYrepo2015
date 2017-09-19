function Fs = Fsam(D);
% adc_data/Fsam - sample frequency of ADC data object.
%   Fsam(D) returns the sample frequency in Hz of adc_data object D.
%
%   See also adc_data/samples, adc_data/samperiod.

Fs = D.Data.TimingCalib.Fsam;











