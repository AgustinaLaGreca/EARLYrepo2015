function Dev = sourceDeviceInfo(Dum, RecInstr); %#ok<INUSL>
% datagrabber/sourceDeviceInfo device info associated with data type
%   Dev = sourceDeviceInfo(datagraber(), RecInstr) returns struct Dev with 
%   information on the source device (and its settings) associated with the
%   recording instructions RecInstr (see experiment/recordingintructions).
%   Each known type of data and/or source device has its own "inside
%   information." Unkown data types produce an error.
%
%   See also datagrabber/allKnowndatatypes, datagrabber/device.

if nargin<2,
    Dev = [];
    return;
end

RecInstr = arginDefaults('RecInstr', []);
if isempty(RecInstr),
    Dev = [];
    return;
end

if ~ismember(RecInstr.DataType, allKnownDataTypes(datagrabber())),
    error(['Unknown data type ''' RecInstr.DataType '''. You may introduce data types by editing Datagrabber/allKnowndatatypes']);
end

Dev.Type = '';
Dev.DataType = RecInstr.DataType;
switch RecInstr.DataType,
    case 'DA-test',
        Dev.Type = 'DAC-output';
        Dev.remark = 'test: DAC output directly connected to ADC input';
    case 'LaserDispl',
        Dev.Type = 'Polytec_OFV5000_Displ';
        Dev.remark = 'Main output of Polytex OFV5000 displacement decoder';
    case 'LaserVelocity',
        Dev.Type = 'Polytec_OFV5000_Velo';
        Dev.remark = 'Main output of Polytex OFV5000 velocity decoder';
    case 'Microphone-1',
        Dev.Type = 'Nexus chan 1';
        Dev.remark = 'B&K microphone connected to chan 1 of Nexus amplifier';
    case 'Microphone-2',
        Dev.Type = 'Nexus chan 2';
        Dev.remark = 'B&K microphone connected to chan 2 of Nexus amplifier';
    case 'Micropipette',
        Dev.Type = '';
        Dev.remark = 'Amplified Round window electrode signal';
    case 'Patch_pipette_1',
        Dev.Type = '';
        Dev.remark = 'Primary output of Chan 1 of Multiclamp 700B';
    case 'Patch_pipette_2',
        Dev.Type = '';
        Dev.remark = 'Primary output of Chan 2 of Multiclamp 700B';
    case 'RW_electrode',
        Dev.Type = '';
        Dev.remark = 'Amplified Round window electrode signal';
    case 'Spikes', 
        Dev.Type = '';
        Dev.remark = 'Spikes, time-stamped by window discriminator or peak picker.';
    case 'Laser reflection', 
        Dev.Type = 'OFV5000';
        Dev.remark = 'Intensity of reflected bundle as number in range [1..256].';
    otherwise,
        error(['Unknown datatype ''' RecInstr.DataType '''']);
end
rs = RecInstr.RecSettings;
if isempty(Dev.Type),
    if isfield(rs, 'Amplifier'), Dev.Type = rs.Amplifier;
    elseif isfield(rs, 'Discriminator'), Dev.Type = rs.Discriminator;
    else, Dev.Type = '???';
    end
end
Dev.RecInstr = RecInstr;
Dev.manualSettings = RecInstr.RecSettings;







