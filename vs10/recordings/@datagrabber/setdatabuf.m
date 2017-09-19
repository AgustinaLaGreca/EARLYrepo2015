function G = setdatabuf(G, Fld, X);
% datagrabber/setdatabuf - set (field of) data buffer of datagrabber object
%   G = setdata(G, 'Foo', X) assigns value X to the Foo field G's DataBuf.
%
%   See also datagrabber/setdatabuf, Grabevents/getdata, Dataset/getdata.

% simply return the Data field
G.DataBuf.(Fld) = X;








