function ET = digichan(DS, Chan);
% Dataset/digichan - digital data of dataset
%    digichan(DS, Chan) returns the events  (an eventdata object)
%    in digital channel Chan. If Chan is a number, say 1, data are obtained 
%    from the RX6_digital_1 field of DS.Data. Alternatively, a full field 
%    name may be specified, e.g., Chan = 'RX8_digital_5'. 
%
%    See also dataset/spiketimes, anachan.

% get the appropriate field of DS.Data
if isnumeric(Chan),
    Chan = ['RX6_digital_' num2str(Chan)]; % default event channel (see help text)
end
ET = DS.Data.(Chan); % eventdata 

