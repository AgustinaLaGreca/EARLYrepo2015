function D = sourcedevice(G, Fld);
% grabbeddata/source - source device of grabbeddata
%   S=sourcedevice(G) returns the info on the source device of G in s
%   Struct S, which contains the fields
%        ID: "fixed info" on the source device, as specified at experiment
%             definition time.
%    Settings: parameter settings polled just prior to the recording (only
%             for devices that can be queried.
%
%   sourcedevice(G, 'Foo') directly returns field Foo of the sourcedevice.
%
%   See also grabbeddata/devicesettings, datagrabber/getsourcedevicesettings.

D = sourcedevice(supplier(G));
if nargin>1,
    D = D.(Fld);
end










