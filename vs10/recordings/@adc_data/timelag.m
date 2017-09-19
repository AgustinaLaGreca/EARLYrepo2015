function t0 = timelag(D);
% adc_data/timelag - time lag of ADC data object.
%   timelag(D) returns the time lag of D, which is the combined AD and DA
%   delay in ms as evaluated prior to the recording by grab_adc/prepare.
%
%   See also adc_data/samples, adc_data/Fsam, adc_data/samperiod,
%   grab_adc/prepare.

t0 = D.Data.TimingCalib.Lag_ms;











