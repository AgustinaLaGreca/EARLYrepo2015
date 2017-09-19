function GD=grabbeddata(SupStatus, Supplier);
% grabbeddata - construct grabbeddata object
%    G = Grabbeddata(SupStatus, S) creates a Grabbeddata object for supplier S.
%
%    Grabbeddata is a superclass of Eventdata, adc_data, comport_data.
%
%    Datagrabber subclasses and Grabbeddata subclasses come in matched
%    pairs, e.g., Grabevents <--> Eventdata.
%
%    See also, Grabevents, Eventdata, dataset/getdata, Datagrabber.

[SupStatus, Supplier] = arginDefaults('SupStatus/Supplier',[]);

GD = CollectInStruct(SupStatus, Supplier);
GD = class(GD, mfilename);


