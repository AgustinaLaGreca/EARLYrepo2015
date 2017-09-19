function G = grab_adc(RecInstr, Hardware, DS, Ax, Name, You);
% grab_adc - grab ADC data
%   G = grab_adc(RecInstr, Hardware, DS, Ax, Name, You) creates a
%   grab_adc Action object G whose task it is to download samples from an
%   ADC channel and to store them in dataset DS.  As long data 
%   collection hasn't finished, the data are kept in G and DS only has 
%   G's address. But DS also gets G's address, thus enabling the final 
%   transfer of the data to DS (at wrapup time), which is initiated by DS.
%   Inputs are
%   RecInstr: detailed recording instructions as provided by
%             Experiment/recordInstructions. The following fields of
%             RecInstr are compulsory:
%          datafieldname: fieldname for storage in dataset DS
%                   Fsam: sample rate in Hz
%           grabInterval: interval in ms between consecutive calls to
%                         the recording device to upload samples. If
%                         grabInterval is an array [dt t0], then dt is the
%                         interval and t0 the waiting time preceding the
%                         first upload call.
%   RecInstr: single element of struct array returned by
%             RecordingInstructions. Contains sample rate, channel
%             specification, grab interval, etc.
%             The datafieldname field of RecInstr is the field name of 
%             the DS.Data subfield where the events wil eventually
%             be stored, e.g., 'RX6_analog_1'. The data will only be transferred
%             to DS at Wrapup time. As long data collection hasn't
%             finished, the data are kept in G and DS only has G's address.
%             Note that multiple grab_adc objects may supply to one and
%             the same dataset, in which case each one needs a unique FN.
%   RecHardware: generic info on hardware settings of the recording device.
%             See recHardware.
%         DS: dataset where data are stored, or its address. DS must have
%             an address, since this is what is stored in G, not DS itself.
%         Ax: address info to store G.
%   The grab_adc is a subclass of both the Action and Datagrabber classes.
%
%   grab_adc(RecInstr, Hardware, DS, FN, figh, 'Foo') is the same as
%   grab_adc(RecInstr, Hardware, DS, FN, GUIaccess(figh, 'Foo')).
%
%   G = grab_adc(RecInstr, Hardware, DS, FN, figh, 'Foo', 'You')
%   makes the finished status of G dependent on the finished status of the
%   action object names 'You' in the same figh GUI (see action/AfterYou).
%   Instead of a single You, a char string of names {'Me' 'You' 'Them'}
%   is also okay, in which case G is ready whenever all of the listed
%   Actions are.
%
%   grab_adc is tightly knot with the adc_data class, the data it is
%   supposed to produce.
%
%   Type 'methods grab_adc -full' to see what can be done with grab_adc
%   objects.
%
%   See also grab_adc/getdata, dataset, Action, Datagrabber,
%   dataset/addsupplier, dataset/getdata, action/AfterYou.

[Name, You] = arginDefaults('Name/You', '');

isVoid=0;
if nargin<1, % void object
    [Fsam, Dev, BufTag, NsamTag, BufSize, iChan, Hardware] = deal([]);  % G's own fields
    [RecInstr, DataName, DataBuf, DS, dt] = deal([]); % superclass constructor args
    isVoid=1;
else,
    Fsam = RecInstr.Fsam; % sample rate in Hz
    iChan = RecInstr.Chan;
    dt = RecInstr.grabInterval; 
    Dev = RecInstr.Dev; % Hardware device from which to grab the events
    DS = address(DS); % access by ref
    % tags in the RPvdS circuit that give access to circuit state
    if isequal(128, iChan),
        BufTag = 'ADC_1'; % the recorded samples, chan 1
        NsamTag = 'Nsam_ADC_1'; % sample count of recorded samples, chan 1
    elseif isequal(129, iChan), 
        BufTag = 'ADC_2';  % ditto chan 2
        NsamTag = 'Nsam_ADC_2'; % ditto chan 2
    else,
        error(['Cannot download samples from ADC channel ' num2str(iChan) '.']);
    end
    BufSize = nan; % will be queried @ prepare time
    % data storage
    DataBuf.Samples = hoard; % the recorded samples themselves are stored dynamically. oneshot adds to it
    DataBuf.TimingCalib = []; % at prepare time: info on calibration of the timing
    DataBuf.FinalStatus = ''; % at wrapup time: a flag indicating status of G at wrapup time
    DataBuf.Supplier = []; % at wrapup time a copy of G itself will here be stored with the data
    DataName = RecInstr.datafieldname;
end

G = CollectInStruct(Fsam, Dev, BufTag, NsamTag, BufSize, iChan, Hardware);
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



