function G = getsourcedevicesettings(G);
% datagrabber/getsourcedevicesettings get settings of source device
%   getsourcedevicesettings(G) checks whether G's source device (the device
%   providing the data) has any settings that can be retrieved over e.g. 
%   a COM port, and stores any settings in the Device property of G.
%   Identification of the device is based on dataType(G).
%
%   For non-queried devices, fixed settings are assumed as described in the
%   code.
%
%   See also datagrabber/datatype.


DT = dataType(G);
s = [];
switch DT,
    case 'DA-test', % stimulus DA directy looped to AD
        s.ScaleFactor = 1;
        s.Unit = 'V';
    case 'LaserVelocity', 
        s = getOFV5000;
        s.b____conversion_________ = '____________';
        s.ScaleFactor = s.VeloRange;
        s.Unit = 'mm/s';
    case 'LaserDispl'
        s = getOFV5000;
        s.b____conversion_________ = '____________';
        s.ScaleFactor = s.DisplRange;
        s.Unit = 'nm';
    case 'Microphone-1', % assuming 0.5" microphone, 1 V/Pa
        s.Microphone = 'B&K 0.5"';
        s.ScaleFactor = 1;
        s.Unit = 'Pa';
    case 'Microphone-2', % assuming 0.25" microphone, 100 mV/Pa
        s.Microphone = 'B&K 0.25"';
        s.ScaleFactor = 0.1;
        s.Unit = 'Pa';
    case {'Micropipette', 'RW_electrode'}
        s.ScaleFactor = nan; % refer to manual settings
        s.Unit = '????';
    case {'Patch_pipette_1', 'Patch_pipette_2'};
        warning('700B communication not yet implemented');
        s.ScaleFactor = nan;
        s.Unit = '????';
    otherwise,
        s.ScaleFactor = nan;
        s.Unit = '????';
end
G.SourceDevice.Settings = s;








