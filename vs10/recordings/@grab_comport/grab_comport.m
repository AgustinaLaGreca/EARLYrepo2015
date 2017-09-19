function G = grab_comport(RecInstr, Hardware, DS, Ax, Name, You);
% grab_comport - grab data from serial port (COM port) of PC
%   G = grab_comport(RecInstr, RecHardware, DS, Ax, Name, You) creates a
%   grab_comport Action object G whose task it is to retrieve information
%   from a device connected to a COM port of the PC and to store it in
%   dataset DS. As long data collection hasn't finished, the data are kept 
%   in G and DS only has G's address. But DS also gets G's address, thus
%   enabling the final "transfer" of the data to DS (at wrapup time), which
%   is initiated by DS.
%   The input arguments are
%   RecInstr: detailed recording instructions as provided by
%             Experiment/recordInstructions. The following fields of
%             RecInstr are compulsory:
%          datafieldname: fieldname for storage in dataset DS
%                   Fsam: sample rate in kHz
%           grabInterval: interval in ms between consecutive calls to
%                         the recording device to upload samples. If
%                         grabInterval is an array [dt t0], then dt is the
%                         interval and t0 the waiting time preceding the
%                         first upload call.
%   RecHardware: generic info on hardware settings of the recording device.
%             See recHardware.
%         DS: dataset where data are stored, or its address. DS must have
%             an address, since this is what is stored in G, not DS itself.
%         Ax: address info to store G.
%   Grab_comport is a subclass of both the Action and Datagrabber classes.
%
%   grab_comport(RecInstr, Hardware, DS, FN, figh, 'Foo') is the same as
%   grab_comport(RecInstr, Hardware, DS, FN, GUIaccess(figh, 'Foo')).
%
%   G = grab_comport(RecInstr, Hardware, DS, FN, figh, 'Foo', 'You')
%   makes the finished status of G dependent on the finished status of the
%   action object names 'You' in the same figh GUI (see action/AfterYou).
%   Instead of a single You, a char string of names {'Me' 'You' 'Them'}
%   is also okay, in which case G is ready whenever all of the listed
%   Actions are.
%
%   grab_comport is tightly knot with the comport_data class, the data it is
%   supposed to produce.
%
%   Type 'methods grab_comport -full' to see what can be done with grab_comport
%   objects.
%
%   See also grab_comport/getdata, dataset, Action, Datagrabber,
%   dataset/addsupplier, dataset/getdata, action/AfterYou.

%     datafieldname: 'PC_COM1'
%               Dev: 'COM1-port'
%           RecType: 'PC_serial'
%           grabber: @grab_comport
%              Chan: 1
%          DataType: 'Laser reflection'
%       RecSettings: []
%      grabInterval: [320 144]
%              Fsam: []

% 'Fsam' 'Dev' 'iChan' 'Hardware' 'DataName'
% 'Data' 'DS' 'datagrabber'    'action'


[Name, You] = arginDefaults('Name/You', '');

isVoid=0;
if nargin<1, % void object
    [Fsam, Dev, Port, PrepCall, QueryCall, CloseCall, TimeCall, Hardware] = deal([]);  % G's own fields
    [RecInstr, DataName, DataBuf, DS, dt] = deal([]); % superclass constructor args
    isVoid=1;
else,
    dt = RecInstr.grabInterval; 
    Fsam = RecInstr.Fsam; % *stimulus* sample rate in Hz
    iChan = RecInstr.Chan;
    Dev = RecInstr.Dev; % Hardware device from which to grab the events
    DS = address(DS); % access by ref
    [PrepCall, QueryCall, CloseCall, TimeCall, Port] = local_getcalls(RecInstr);
    % data storage
    DataBuf.Samples = hoard; % the recorded samples themselves are stored dynamically. oneshot adds to it
    DataBuf.Times = hoard; % the times at which the samples are retrieved. same procedure as teh samples.
    DataName = RecInstr.datafieldname;
end

G = CollectInStruct(Fsam, Dev, Port, PrepCall, QueryCall, CloseCall, TimeCall, Hardware);
G = class(G,mfilename, datagrabber(RecInstr, DataName, DataBuf, DS), action('initialized',dt));
% now the object exists, but we're not done yet
if ~isVoid, % upload G and subscribe G as a supplier of dataset DS
    if isempty(Name), ADR={Ax}; else, ADR = {Ax Name}; end
    [G, adG] = upload(G, ADR{:});
    % subscribe as a supplier of dataset DS.
    addsupplier(download(DS), adG, DataName); % download(DS) needed to reach dataset method addsupplier ..
    % if a You is specified, define finish dependency
    if ~isempty(You), G=afteryou(G, Ax, You); end
end;

%=======================


function [PrepCall, QueryCall, CloseCall, TimeCall, Port] = local_getcalls(RI);
% ugly .. recordingsinstructions should be made more flexible
Port = ['COM' num2str(RI.Chan)]; % COM1 etc
switch RI.DataType, 
    case 'Laser reflection',
        PrepCall = {@getOFV5000 Port 'SignalLevel' 1}; % 1: keep conn open
        QueryCall = {@getOFV5000 Port 'SignalLevel' 1};
        TimeCall = {@seqplaystatus 1 'isam_abs'};
        CloseCall = {@getOFV5000 'close'};
    otherwise, error(['Unknown serial port data ''' RI.DataType, '''.']);
end





