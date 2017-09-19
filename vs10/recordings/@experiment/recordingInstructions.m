function [RI, CI] = recordingInstructions(E, StimParam, RecParam);
% recordingInstructions - assemble instructions for a recording
%   [RI, CI]=recordingInstructions(Exp, StimParam, RecParam) returns
%   recording instructions RI and circuit instructions CI to be used
%   by PlayRecStop at recording time.
%   Inputs:
%           Exp: experiment definition containing main recording parameters
%     StimParam: stimulus parameters from stimulus GUI
%      RecParam: additional recording parameters from dashboard GUI.
%      
%   RecordingInstructions is typically called by Dashboard when the Play
%   Record button are pushed. 
%
%   See also StimGUI, recGUIpanel, Experiment.

RX6_DADig_timing = 0; % default: DA->Dig timing calibration of RX6 not needed
RX6_DAAD_timing = 0; % default: DA->AD timing calibration of RX6 not needed

RI = [];
Erec = E.Recording.General;

% get list of all active recording channels specified to be used for the
% current recording.
[RecSrc, SrcNames] = recordingsources(E,1,RecParam);

% ========event timing over Dig0 ("A1") of RX6
Chan = 'RX6_digital_1';
if ismember(Chan, SrcNames),
    src = RecSrc.(Chan);
    ri = struct('datafieldname', Chan);
    ri.Dev = 'RX6';
    ri.RecType = 'RX6_digital';
    ri.grabber = @grabevents;
    ri.Chan = 0;
    ri.DataType = src.DataType;
    ri.RecSettings.Discriminator = src.Discriminator;
    ri.grabInterval = [133 400]; % ms grab interval; wait a while before starting grabbing
    RI = [RI ri];
    RX6_DADig_timing = 1; % we need DA->Dig timing calibration of RX6
end

% ========analog recording over chan-1
Chan = 'RX6_analog_1';
if ismember(Chan, SrcNames),
    src = RecSrc.(Chan);
    ri = struct('datafieldname', Chan);
    ri.Dev = 'RX6';
    ri.RecType = 'RX6_analog';
    ri.grabber = @grab_adc;
    ri.Chan = 128;
    ri.DataType = src.DataType;
    ri.RecSettings.Amplifier = src.Amp;
    ri.RecSettings.Gain = src.Gain;
    ri.RecSettings.GainUnit = src.GainUnit;
    ri.RecSettings.lowCutoff = src.LowFreq;
    ri.RecSettings.highCutoff = src.HighFreq;
    ri.grabInterval = [300 1500]; % ms grab interval; start delay of 1500 ms (to give the other actions some air)
    RI = [RI ri];
end

% ========analog recording over chan-2
Chan = 'RX6_analog_2';
if ismember(Chan, SrcNames),
    src = RecSrc.(Chan);
    ri = struct('datafieldname', Chan);
    ri.Dev = 'RX6';
    ri.RecType = 'RX6_analog';
    ri.grabber = @grab_adc;
    ri.Chan = 129;
    ri.DataType = src.DataType;
    ri.RecSettings.Amplifier = src.Amp;
    ri.RecSettings.Gain = src.Gain;
    ri.RecSettings.GainUnit = src.GainUnit;
    ri.RecSettings.lowCutoff = src.LowFreq;
    ri.RecSettings.highCutoff = src.HighFreq;
    ri.grabInterval = [320 1650]; % ms grab interval; start delay of 1650 ms (to give the other actions some air)
    RI = [RI ri];
end

% ========analog recording over chan-2
Chan = 'PC_COM1';
if ismember(Chan, SrcNames),
    src = RecSrc.(Chan);
    ri = struct('datafieldname', Chan);
    ri.Dev = 'COM1-port';
    ri.RecType = 'PC_serial';
    ri.grabber = @grab_comport;
    ri.Chan = 1;
    ri.DataType = src.DataType;
    ri.RecSettings = []; % to be queried
    ri.grabInterval = [120 120]; % ms grab interval; start after 120 ms.
    RI = [RI ri];
end

% RX6 circuit info
Fsam = StimParam.Fsam; % RX6 stimulus sample rate in Hz
if ~isempty(RI), 
    [RI.Fsam] = deal(Fsam); % this is the stimulus sample rate, which can serve as master clock
end % conditional because it is not certain if any recording is requested
if Fsam<minFsam(E), error('Stimulus sample rate is smaller than sample rate dictated by experiment specs.'); end
if isequal('Left', AudioChannelsUsed(E)),
    %changed by Abel on 3/05/2016
    %Using monaural (Left play/Left rec) renders no spikes and no errors
    %RX6seqplay-trig-2ADC-DA1only.rcx doesn't seem to work
    %TODO: Ask Marcel, why -DA1only in this case? (maybe because the time
    %calibration is done on the left channel?)
    RX6_circuit = getFieldOrDefault(StimParam, 'RX6_circuit', 'RX6seqplay-trig-2ADC');
    %Original value:
    %RX6_circuit = getFieldOrDefault(StimParam, 'RX6_circuit', 'RX6seqplay-trig-2ADC-DA1only');
else,
    RX6_circuit = getFieldOrDefault(StimParam, 'RX6_circuit', 'RX6seqplay-trig-2ADC');
end
N_DSP = sys3setup('DSPcount_RX6'); % # DSPs in RX6
if isequal(2, N_DSP), % default is 5 DSPs. Use special, "cheap" circuit when only 2 DSPs are available
    RX6_circuit = [RX6_circuit '-2DSP']; 
end
RX6_Fsam = Fsam;

% RP2 circuit info (not implemented yet)
RP2_circuit = '';
RP2_Fsam = nan;
createdby = mfilename;

CI = CollectInStruct(RX6_circuit, RX6_Fsam, RX6_DAAD_timing, RX6_DADig_timing, '-', ...
    RP2_circuit, RP2_Fsam, createdby);




