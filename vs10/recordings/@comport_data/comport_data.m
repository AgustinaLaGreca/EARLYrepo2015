function D = comport_data(Data, SupStatus, Supplier);
% comport_data - construct comport_data object.
%    D=comport_data(Data, SupStatus, Supplier) is an comport_data object whose
%    data are in Data and whose supplier is Supplier. SupStatus is the
%    supplier status.
%    Not intended for direct use; only to be called by grab_comport/getdata.
%    comport_data is a subclass of Grabbeddata.
%
%    See also grabevents/getdata, Grabbeddata, dataset/getdata, 
%    comport_data/supplierstatus.


if nargin<1,
    [Data, SupStatus, Supplier] = deal([]);
end
D = CollectInStruct(Data);
D = class(D,mfilename, grabbeddata(SupStatus, Supplier));





