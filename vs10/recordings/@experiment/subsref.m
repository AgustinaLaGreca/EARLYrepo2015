function Y = subsref(E, S);
% experiment/SUBSREF - SUBSREF for Experiment objects
%
%   E(5) is the 5th element of Experiment array E.
%
%   E{..} is not a valid syntax.
%
%   E.Foo returns property Foo of E. These properties directly reflect the 
%   specs provided by the user during the use of experiment/edit. Here is a 
%   list of valid properties (unlike MATLAB fieldnames, they are case INsensitive):
%  --HELP--
%              Help: isplay this help text.
%   --ID properties--
%              Name: name of experiment
%      Experimenter: name of experimenter (see experiment/isknownexperimenter)
%              Type: experiment type
%           Started: time of creation in datestr format
%          Computer: computer at which the experiment was defined (see compuname)
%          Location: ditto location (see where)
%
%   --Audio properties--
%          MinStimFreq: min stimulus freq [Hz]
%          MaxStimFreq: max stimulus freq [Hz]
%    AudioChannelsUsed: DA channels used for sound stimulation
%          Attenuators: attenuators used
%    PreferredNumAtten: preferred minimum numerical attenuation
%      CalibrationType: calibration type
%               Probes: probe tubes used
%
%   --Recording properties--
%      RecordingSide: which side of the animal do we record from
%          EventType: what type of events are recorded
%      Discriminator: which spike discriminator is used
%         AD1_signal: what type of signal is recorded over chan 1 of A/D
%      AD1_Amplifier: which conditioning amplifier is used for AD(1)
%           AD1_Gain: what is its gain or sensitivity
%       AD1_GainUnit: in what units is the gain/sensitivity specfied
%      AD1_lowCutoff: low cutoff [Hz] of filter
%     AD1_highCutoff: high cutoff [Hz] of filter
%         AD2_signal: what type of signal is recorded over chan 2 of A/D
%      AD2_Amplifier: which conditioning amplifier is used for AD(1)
%           AD2_Gain: what is its gain or sensitivity
%       AD2_GainUnit: in what units is the gain/sensitivity specfied
%      AD2_lowCutoff: low cutoff [Hz] of filter
%     AD2_highCutoff: high cutoff [Hz] of filter
%
%   --Preferences--
%         maxNcond: max number of condition in a recording
%    ITDconvention: what it means for ITDs to be positive
%  StoreComplexWaveforms: do or don't store complex analytic waveforms
% 
%   --Sys3setup--
%          
%   --Devices--
%        AudioDevice: device used for sound D/A
%      AudioMaxAbsDA: max abs value of samples sent to audio device
%      ElectraDevice: device used for electrical stimulation
%    ElectraMaxAbsDA: max abs value of samples sent to electra device
%    ElectraMinFsam_kHz: min sample rate [kHz] for elactrical stimulation
%    RecordingDevice: device used for recordings and/or event timing
%
%   Subsref is used for all these properties to avoid a proliferation of 
%   trivial "get" methods. More subtle properties have their own method,
%   e.g., Experiment/minFsam.
%
%   Staggered subsref usage as in E.Foo(4).Faa is handled recursively.
%
%   See also dataset, "methods dataset".


if length(S)>1, % use recursion from the left
    y = subsref(E,S(1));
    Y = subsref(y,S(2:end));
    return;
end
%----single-level subsref from here (i.e., S is scalar) ---
switch S.type,
    case '()', % array element
        try, Y = builtin('subsref',E,S);
        catch, error(lasterr);
        end
    case '{}',
        error('Experiment objects do not allow subscripted reference using braces {}.');
    case '.',
        Field = S.subs; % field name
        switch lower(Field),
            %--Help--
            case 'help', help('experiment/subsref');
            %--ID properties--
            case 'id', Y = E.ID;
            case 'name',  Y = E.ID.Name;
            case 'experimenter', Y = E.ID.Experimenter;
            case 'type', Y = E.ID.Type;
            case 'started', Y = E.ID.Started;
            case 'computer', Y = E.ID.Computer;
            case 'location', Y = E.ID.Location;

            %--Audio properties--
            case 'audio', Y = E.Audio;
            case 'minstimfreq', Y = E.Audio.MinStimFreq;
            case 'maxstimfreq', Y = E.Audio.MaxStimFreq;
            case 'audiochannelsused',  Y = E.Audio.DAchannelsUsed;
            case 'attenuators',  Y = E.Audio.Attenuators;
            case 'preferrednumatten',  Y = E.Audio.PreferredNumAtten;
            case 'calibrationtype',  Y = E.Audio.CalibrationType;
            case 'probes',  Y = E.Audio.Probes;
                
            %-- Recording props--
            case 'recording', Y = E.Recording;
            case 'recdevice', Y = E.Recording.General.RecordingSide;
            case 'recminFsam', Y = E.Recording.General.RecordingSide;
            case 'recminfsamunit', Y = E.Recording.General.RecordingSide;
            case {'recordingside', 'recside'}, Y = E.Recording.General.Side;
            case 'source', Y = E.Recording.Source;

            %-- Sys3setup props--
            case 'sys3setup', Y = E.Sys3setup;

            %--Preferences--
            case 'maxncond', Y = E.Preferences.MaxNcond;
            case 'itdconvention', Y = E.Preferences.ITDconvention;
            case 'storecomplexwaveforms', Y = E.Preferences.StoreComplexWaveforms;

            %--Devices--
            case 'audiodevice', Y = E.Audio.Device;
            case 'audiomaxabsda', Y = E.Audio.MaxAbsDA;
            case 'electradevice', Y = E.Electra.Device;
            case 'electramaxabsda', Y = E.Electra.MaxAbsDA;
            case 'electraminfsam_khz', Y = E.Electra.MinFsam_kHz;
            case 'recordingdevice', Y = E.Recording.General.Device;
            otherwise,
                error(['Invalid Experiment property ''' Field '''. See Experiment/subsref for a list of valid properties.'])
        end
end % switch/case


