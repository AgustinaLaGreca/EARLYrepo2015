function D = dataType(G, varargout);
% datagrabber/dataType - data type of datagrabber
%   dataType(G) returns a char string specifying the type of data stored in
%   it.
%
%   See also grabbeddata.

RI = recordinstructions(G);
D = RI.DataType;









