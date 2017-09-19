function D = eventdata(Data, SupStatus, Supplier);
% eventdata - construct eventdata object.
%    D=Eventdata(Data, SupStatus, Supplier) is an Eventdata object whose
%    data are in Data and whose supplier is Supplier. SupStatus is the
%    supplier status.
%    Not intended for direct use; only to be called by grabevents/getdata.
%    Eventdata is a subclass of Grabbeddata.
%
%    See also grabevents/getdata, Grabbeddata, dataset/getdata, 
%    eventdata/supplierstatus.


if nargin<1,
    [Data, SupStatus, Supplier] = deal([]);
end
D = CollectInStruct(Data);
D = class(D,mfilename, grabbeddata(SupStatus, Supplier));





