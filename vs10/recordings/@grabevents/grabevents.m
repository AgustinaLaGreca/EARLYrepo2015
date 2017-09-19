function G = grabevents(RecInstr, Hardware, DS, Ax, Name, You);
% grabevents - grab timed events
%   G = grabevents(RecInstr, Hardware, DS, Ax) creates a
%   grabevents Action object G whose task it is to grab timed events from a
%   hardware timer and to store them in dataset DS.  As long data 
%   collection hasn't finished, the data are kept in G and DS only has 
%   G's address. But DS also gets G's address, thus enabling the final 
%   transfer of the data to DS (at wrapup time), which is initiated by DS.
%   Inputs are
%   RecInstr: detailed recording instructions as provided by
%             Experiment/recordInstructions. The following fields of
%             RecInstr are compulsory:
%          datafieldname: fieldname for storage in dataset DS
%                   Fsam: circuit sample rate in Hz
%           grabInterval: interval in ms between consecutive calls to
%                         the recording device to upload samples. If
%                         grabInterval is an array [dt t0], then dt is the
%                         interval and t0 the waiting time preceding the
%                         first upload call.
%   Hardware: recHardware object containing low-level params for recording 
%             device.
%         DS: dataset where data are stored, or its address. DS must have 
%             an address, since this is what is stored in G, not DS itself.
%         Ax: address info to store G itself. 
%   Grabevents is a subclass of both the Action and Datagrabber classes.
%
%   grabevents(RI, Hw, DS, FN, figh, 'Foo') is the same as
%   grabevents(RI, Hw, DS, FN, GUIaccess(figh, 'Foo')), i.e., the
%   grab-events object is stored as Foo property of a figure's GUIdata (see
%   getGUIdata, GUIaccess, Address).
%
%   G = grabevents(R, H, DS, FN, figh, 'Foo', 'You') makes the
%   finished status of G dependent on the finished status of the action
%   object names 'You' in the same figh GUI (see action/AfterYou).
%   Instead of a single You, a char string of names {'Me' 'You' 'Them'}
%   is also okay.
%
%   Grabevents is tightly knot with the EventData class, the data it is
%   supposed to produce.
%
%   Type 'methods grabevents -full' to see what can be done with Grabevents 
%   objects.
%
%   See also grabevents/getdata, dataset, Action, Datagrabber,
%   dataset/addsupplier, dataset/getdata, action/AfterYou.

[Name, You] = arginDefaults('Name/You', '');

isVoid=0;
if nargin<1, % void object
    [Fsam, Dev, StampBit, Hardware] = deal([]); % G's own fields
    [RecInstr, DataName, DataBuf, DS, dt] = deal([]); % superclass constructor args
    isVoid=1;
else,
    DS = address(DS); % access by ref
    Dev = Hardware.DAC; % TDT device from which to grab the events
    StampBit = Hardware.SpikeTimeStampBit;
    % extract basic info from recording instructions
    DataName = RecInstr.datafieldname;
    Fsam = RecInstr.Fsam;
    dt = RecInstr.grabInterval;
    % data storage: specific to event type data
    DataBuf.EventTimes = hoard([],inf); % dynamic storage of events. Oneshot gets data. inf: no streaming to disk
    DataBuf.TimingCalib = []; % at prepare time: info on calibration of the timing
    DataBuf.FinalStatus = ''; % at wrapup time: a flag indicating status of G at wrapup time
    DataBuf.Supplier = []; % at wrapup time a copy of G itself will here be stored with the data
end

G = CollectInStruct(Fsam, Dev, StampBit, Hardware);
G = class(G,mfilename, datagrabber(RecInstr, DataName, DataBuf, DS), action('initialized',dt));
% now the object exists, but we're not done yet
if ~isVoid, % upload G and subscribe G as a supplier of dataset DS
    if isempty(Name), ADR={Ax}; else, ADR = {Ax Name}; end
    [G, adG] = upload(G, ADR{:}); 
    % subscribe as a supplier of dataset DS.
    addsupplier(download(DS), adG, DataName); % download(DS) needed to reach dataset method addsupplier ..
    % if a You is specified, define a finish dependency
    if ~isempty(You), G=afteryou(G, Ax, You); end
end;





