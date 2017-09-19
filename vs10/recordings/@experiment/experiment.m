function E = experiment(varargin);
% experiment - construct experiment object
%    Experiment() returns a void Experiment object, which can be filled
%    using Experiment/Edit.
%
%    Experiment(S), where S is a struct,  converts S to a an Experiment object.
%
%    Experiment(E), where E is an Experiment object, returns E.
%
%    To register an experiment (create a directory, store the status, etc)
%    use Experiment/checkin.
%    To retrieve a named experiment, see Experiment/load.

Version = 1; % initial version
ID = local_ID(Version);
Audio = local_Audio(Version);
Electra = local_Electra(Version);
Recording = struct('General', [], 'Source', []);
Preferences = local_Preferences(Version);
Sys3setup = local_Sys3setup(Version);
Status = local_Status(Version);

Template = CollectInStruct(ID, Audio, Electra, Recording, Sys3setup, ...
    Preferences, Status, Version);

if nargin<1, % Experiment(): void Experiment object
    E = Template;
elseif isstruct(varargin{1}), % Experiment(S): convert struct S to Experiment object
    E = union(Template, varargin{1}); % incorporate any new fields added since previous version upgrade
elseif isa(varargin{1},mfilename), % Experiment(E): tautology
    E = varargin{1};
    return;
end

E = class(E,mfilename);

%=============================
function V = local_ID(Version);
switch Version
    case 1,
        V = VoidStruct('Name', 'Experimenter', 'Type', 'Started', ...
            'Computer', 'Location');
    otherwise, error(['Invalid version ' num2str(Version)]);
end
%=============================
function A = local_Audio(Version);
switch Version
    case 1,
        A = VoidStruct('Device', 'MaxAbsDA', 'DAchannelsUsed', ...
            'HeadphoneBuffer', 'Speakers', 'Probes', ...
            'Attenuators', 'PreferredNumAtten', ...
            'MinFsam_kHz', 'MinStimFreq', 'MaxStimFreq', ...
            'CalibrationType');
    otherwise, error(['Invalid version ' num2str(Version)]);
end

%=============================
function V = local_Electra(Version);
switch Version
    case 1,
        V = VoidStruct('Device', 'DA1_target', 'DA2_target', 'MaxAbsDA', 'MinFsam_kHz');
    otherwise, error(['Invalid version ' num2str(Version)]);
end
%=============================
function V = local_Preferences(Version);
switch Version
    case 1,
        V = VoidStruct('MaxNcond', 'ITDconvention', 'StoreComplexWaveforms', 'CheckID','CheckFullDsInfo');
    otherwise, error(['Invalid version ' num2str(Version)]);
end

%=============================
function V = local_Status(Version);
switch Version
    case 1,
        V = VoidStruct('Modified', 'StatusModified', 'Finished',  ...
            'State', 'Ndataset', 'IDlastSaved', ...
            'iPen', 'PenDepth', 'iCell', 'iRecOfCell', ...
            'AllSaved', 'iCalibMeasured');
    otherwise, error(['Invalid version ' num2str(Version)]);
end

function V = local_Sys3setup(Version);
try,
    switch Version
        case 1,
            V = sys3setup;
        otherwise, error(['Invalid version ' num2str(Version)]);
    end
catch
    V = [];
end




