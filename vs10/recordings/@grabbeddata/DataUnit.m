function U = DataUnit(G);
% grabbeddata/DataUnit - physical unit of grabbed data
%   DataUnit(G) returns a string describing the physical units describe
%   of the recording in G. If no information is know about the physical 
%   units then Unit = 'V', indicating the Voltage of the DAC.
%
%   See also grabbeddata/sourcedevice, grabbeddata/conversionfactor.

[dum, U] = conversionfactor(G);









