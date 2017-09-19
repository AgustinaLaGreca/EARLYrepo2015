function D = adc_data(Data, SupStatus, Supplier);
% adc_data - construct adc_data object.
%    D=adc_data(Data, SupStatus, Supplier) is an adc_data object whose
%    data are in Data and whose supplier is Supplier. SupStatus is the
%    supplier status.
%    Not intended for direct use; only to be called by grab_adc/getdata.
%    adc_data is a subclass of Grabbeddata.
%
%    See also grab_adc/getdata, Grabbeddata, dataset/getdata, 
%    adc_data/supplierstatus.

if nargin<1,
    [Data, SupStatus, Supplier] = deal([]);
end
D = CollectInStruct(Data);
D = class(D,mfilename, grabbeddata(SupStatus, Supplier));
