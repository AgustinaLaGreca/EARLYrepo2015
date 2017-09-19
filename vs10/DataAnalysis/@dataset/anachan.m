function AD = anachan(DS, Chan);
% Dataset/anachan - analog data of dataset
%    anachan(DS, Chan) returns the AD recording (an adc_data object)
%    in analog channel Chan. If Chan is a number, say 2, data are obtained 
%    from the RX6_analog_2 field of DS.Data. Alternatively, a full field 
%    name may be specified, e.g., Chan = 'RX8_analog_13'. 
%
%    See also adc_data/Samples, digichan.

% get the appropriate field of DS.Data
if isnumeric(Chan),
    Chan = ['RX6_analog_' num2str(Chan)]; % default AD channel (see help text)
end
AD = DS.Data.(Chan); % analog data 

