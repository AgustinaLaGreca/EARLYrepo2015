function H = rechardware(Exp, Circuit)
% recHardware - hardware specs for recordings
%    rechardware(Exp, Circuit) returns a object containing specs of the
%    hardware used for recordings. Exp is the Experiment definition, which
%    includes the RecordingDevice property (see Experiment/subsref).
%    Circuit is the name of TDT circuit used for simultaneous D/A and A/D.
%    Currently, only the RX6 and Seqplay** circuits are supported.

if nargin<2,
    Circuit = FromSetupFile('seqplay', 'DefaultCircuitPostfix', '-default', '');
end

if nargin<1,
    H = VoidStruct('DAC', 'SeqplayEnabled', 'SeqplayPostFix', 'SpikeTimeStampBit', 'CalibTimeStampBit', ...
        'ResetDAC_trigger', 'CalibStamping_trigger', 'ResetStamping_trigger', 'CalibADC_trigger', 'ResetADC_trigger', ...
        'AttenType', 'MaxAtten', 'PrefMinNumAtten');
else,
    if isvoid(Exp), % use default RX6/PA5/10 dB
        clear Exp;
        Exp.RecordingDevice = 'RX6';
        Exp.Attenuators = 'PA5';
        Exp.PreferredNumAtten = 10; % dB
    end
    switch upper(Exp.RecordingDevice),
        case 'RX6',
            H.DAC = 'RX6';
            H.SeqplayEnabled =  1;
            H.SeqplayPostFix = strrep(Circuit, 'RX6seqPlay', ''); % only the postfix
            H.SpikeTimeStampBit = 0;
            H.CalibTimeStampBit = 3;
            H.ResetDAC_trigger = @()sys3trig(1,'RX6');
            H.CalibStamping_trigger = @()sys3trig(2,'RX6');
            H.ResetStamping_trigger = @()sys3trig(3,'RX6');
            H.CalibADC_trigger = @()sys3trig(4,'RX6');
            H.ResetADC_trigger = @()sys3trig(5,'RX6');
            H.AttenType = Exp.Attenuators;
            H.MaxAtten = 120;
            H.PrefMinNumAtten = Exp.PreferredNumAtten;
        otherwise,
            error(['Recording hardware device ''' Exp.RecordingDevice ''' not yet supported.'])'
    end
end
H = class(H, mfilename);









