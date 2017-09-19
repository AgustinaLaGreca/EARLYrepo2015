function DG=datagrabber(RecInstr, DataName, DataBuf, DS);
% datagrabber - construct datagrabber object
%    Datagrabber(RecInstr, DataName, DataBuf, DS) creates a Datagrabber object 
%    which grabs data according to recording instructions RecInstr and 
%    stashes them in in data buufer DataBuf. The final destination 
%    of the data is dataset DS, where the Data will be stored under the 
%    name DataName.
%
%    Datagrabber is a superclass to which Action objects like
%    Grabevents, Grab_adc, and grab_comport belong. Datagrabber objects 
%    have Getdata, Clear, and Stop methods, but it is an error not to 
%    overload these methods for any concrete subclass of Datagrabber.
%
%    See also Grabevents, dataset/getdata, Grabbeddata, datagrabber/prepare.

[RecInstr, DataName, DataBuf, DS] = arginDefaults('RecInstr/DataName/DataBuf/DS', []);

if isempty(RecInstr), 
    [SourceDevice, StopStatus] = deal([]);
else, % test if this device is know, andget some additional info on it.
    SourceDevice.ID = sourceDeviceInfo(datagrabber(), RecInstr); % what the device is
    SourceDevice.Settings = []; % device settings may be polled or changed at prepare time.
    StopStatus = [];
end
% 

DG = class(CollectInStruct(RecInstr, SourceDevice, StopStatus, DataName, DataBuf, DS), mfilename);
superiorto('action'); % this forces the overloading of clear; see datagrabber/clear


