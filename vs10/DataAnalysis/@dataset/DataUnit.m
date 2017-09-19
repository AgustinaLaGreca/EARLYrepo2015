function U = DataUnit(DS, Chan);
% dataset/DataUnit - physical unit of recorded data
%    DataUnit(DS, Chan) returns the physical units of the data in channel
%    Chan of dataset DS. If no information is know about the physical 
%   units then Unit = 'V', indicating the Voltage of the DAC.
%
%    See also grabbeddata/DataUnit, dataset/anachan.

U = DataUnit(anachan(DS, Chan));