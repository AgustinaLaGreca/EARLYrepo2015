function addsupplier(D, Address, Trade);
% Dataset/addsupplier - add a data supplier to a dataset object
%    addsupplier(D, adSup, 'Foo') subscribes a supplier with address adSup
%    to dataset D. A field Foo in the Data field of DS is created, and
%    assigned the value adSup. The suplier must be a Datagrabber object.
%
%    Addsupplier is typically invoked by the constructors of datagrabbers
%    like Grabevents. The suppliers "subscribe themselves" to a dataset D.
%    Whereas DataGrabbers like GrabEvents independently dowload data from
%    the hardware, the initiative of gathering these data into a single
%    Dataset object is with the Dataset object (see dataset/getdata).
%
%   See Dataset, Grabevents, grabevents/getdata, dataset/getdata,
%   Datagrabber.

D = download(D);
% check if candidate supplier is a recognized Datagrabber
if ~isa(download(Address),'datagrabber'),
    error('Dataset supplier must be Datagrabber objects.');
end
% create Data subfield - must be unique
if isfield(D.Data, Trade),
    error(['Data subfield named ''' Trade ''' already exists in dataset.']);
end
D.Data(1).(Trade) = Address; % 
upload(D);


