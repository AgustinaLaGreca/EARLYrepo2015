function N = Nsam(D);
% adc_data/Nsam - number of samples of ADC data object.
%   Nsam(D) returns the total number of samples of adc_data object D.
%
%   See also adc_data/samples.

N = sum(D.Data.Samples.NsamSaved);











