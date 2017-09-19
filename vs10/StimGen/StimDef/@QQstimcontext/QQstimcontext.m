function S = stimcontext(DAchan, minFreq, maxFreq, maxNcond, recSide, ITDconvention);
% stimcontext - constructor for stimcontext objects
%    S = stimcontext(DAchan, minFreq, maxFreq, maxNcond, recSide, ITDconvention);
%    StimContext contains stimulus context, i.e., the implicit stimulus
%    parameters that co-determine the stimulus waveform in addition to the
%    user-pecified parameters.
%    ---UNDER CONSTRUCTION---
%
%    ITDconvention = IpsiLeading|ContraLeading|LeftLeading|RightLeading
%                    indicates what positive ITD values mean.

% special cases: struct input or object input
if nargin==1 && isstruct(DAchan),
    S = DAchan;
    S = class(S, mfilename);
    return;
elseif nargin==1 && isa(DAchan, mfilename),
    S = DAchan;
    return;
elseif nargin<1, % void object with correct fields
    [DAchan, minFreq, maxFreq, minADCrate, maxNcond, recSide, ITDconvention, Hardware, Preferences] = deal([]);
    S = CollectInStruct(DAchan, minFreq, maxFreq, minADCrate, maxNcond, recSide, ITDconvention, Hardware, Preferences);
    S = class(S, mfilename);
    return;
end

%------regular call from here: stimcontext is fully specified ------------
[DAchan, Mess] = keywordMatch(DAchan,{'Left' 'Right' 'Both'}, 'DAchan input argument');
error(Mess);
if nargin<2, minFreq=50; end % Hz min stim freq
if nargin<3, maxFreq=45e3; end; % Hz max stim freq
if nargin<4, maxNcond = 400; end; % max # conditions 
if nargin<5, recSide = 'Left'; end; % recording side
if nargin<6, ITDconvention = 'IpsiLeading'; end; % ITD sign convention: what does ITD>0 mean?
% minimum sample rate dictated by recordings over ADC channel (see experiment) 
minADCrate = 0;%minFsam(current(experiment));
[recSide, Mess] = keywordMatch(recSide, {'Left' 'Right'}, 'RecSide input argument');
error(Mess);
[ITDconvention, Mess] = keywordMatch(ITDconvention, ...
    {'IpsiLeading' 'ContraLeading' 'LeftLeading' 'RightLeading'}, 'ITDconvention input argument');
error(Mess);

error(numericTest(maxNcond,'posint','maxNcond arg is '));

Hardware.DAC = 'RX6';
Hardware.SeqplayEnabled = 1;
%Hardware.SeqplayPostFix = '-triggering'; % default postfix of Seqplay rcx circuit.
Hardware.SeqplayPostFix = '-trig-2ADC'; % default postfix of Seqplay rcx circuit.
Hardware.SpikeTimeStampBit = 0; % bit number for digital input of timing the action potentials (spikes)
Hardware.CalibTimeStampBit = 3; % bit number for digital input of timing calibration pulse train (see TimingCalibrate)
% triggers
Hardware.ResetDAC_trigger = @()sys3trig(1,'RX6'); % software trigger to reset D/A conversion 
Hardware.CalibStamping_trigger = @()sys3trig(2,'RX6'); % software trigger to calibrate DAC vs digital in (for event grabbing)
Hardware.ResetStamping_trigger = @()sys3trig(3,'RX6'); % software trigger to reset time stamp buffers for event grabbing
Hardware.CalibADC_trigger = @()sys3trig(4,'RX6'); % software trigger to calibrate DAC vs ADC in (for analog recordings)
Hardware.ResetADC_trigger = @()sys3trig(5,'RX6'); % software trigger to reset ADC buffers for analog recordings

% Attenuators, if present
if ismember('PA5_1', sys3devicelist) && ismember('PA5_2', sys3devicelist),
    Hardware.AttenType = 'PA5';
    Hardware.MaxAtten = 120; % dB maximum analog attenuation
else,
    Hardware.AttenType = 'NONE';
    Hardware.MaxAtten = 0; % dB maximum analog attenuation
end
Hardware.PrefMinNumAtten = 10; % dB preferred minimum numerical attenuation
Preferences.StoreComplexWaveforms = 0;
S = CollectInStruct(DAchan, minFreq, maxFreq, minADCrate, maxNcond, recSide, ITDconvention, Hardware, Preferences);
S = class(S, mfilename);



    
    


