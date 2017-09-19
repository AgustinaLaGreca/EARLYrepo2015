function D = dataset(Experiment, Stim, Rec, varargin);
% dataset - dataset constructor
%   DS = dataset(Experiment, Stim, Rec, Ax) initializes a dataset object.
%   Experiment is an experiment object specifying to which experiment DS 
%   belongs. Stim is the stimulus used to collect DS. Rec is a struct
%   containing recording instructions as returned by RecordingInstructions.
%   Datasets are dynamic objects. The address of DS ("upload destination") 
%   is specified by access struct Ax (see Dynamic). 
%
%   Filling DS with data is done by action ojects of class Supplier.
%
%   DS = dataset(Experiment, Stim, Rec, figh, 'Foo') is the same as
%   DS = dataset(Experiment, Stim, Rec, DataType, GUIaccess(figh,'Foo')).
%
%   Type 'methods dataset -full' to see what can be done with datasets. 
%
%   As of oktober 2016 a new field is added to the Stim struct: Stutter
%   If this field is set to 'On' the first condition is played twice, if it
%   is 'Off' the stimulus is played as usual. (By Adriaan Lambrechts)
%
%   See also Dynamic, Action, Supplier, GUIaccess, Experiment, 
%   RecordingInstructions.

isVoid = 0;
if nargin<1,
    isVoid = 1;
elseif nargin==1 && ischar(Experiment),
    D = dataset();
    D.ID = Experiment;
    return;
elseif nargin<4,
    error('Non void dataset object must have an upload destination.');
end

if isVoid,
    [ID, Stim, Rec, Data, Status] = deal([]);
else,
    ID = local_ID(Experiment);
%     Stim.Waveform = strip(Stim.Waveform); % remove samples from waveforms to save space
    Data = struct([]); % to be filled by the suppliers
    Status = 'hopingfordata';  % no data yet
end

Stopped = [];
D = CollectInStruct(ID, Status, Stim, Rec, Data, Stopped);
D = class(D,mfilename, dynamic);
if ~isVoid, D=upload(D,varargin{:}); end

%==============================
function ID=local_ID(Exp);
% (under development)
% stores info identifying the dataset within the current experiment, etc
ID.Experiment = Exp;
ID.Computer = CompuName;
ID.Location = Where;
ID.Experimenter = username;
ID.created = datestr(now);
ID.modified = datestr(now);
ID.iDataset = nan;
ID.iCell = nan;
ID.iRecOfCell = nan;
ID.iPen = nan;
ID.PenDepth = nan;
rand('twister',sum(100*clock)); ID.uniqueID = RandomInt(1e15);
ID.contributors = {fullfile(mfilename, mfilename('class'))};





