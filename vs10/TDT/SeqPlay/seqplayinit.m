function [y, RCLD] = SeqplayInit(Fsam, Dev, Postfix, Recycle);
% SeqplayInit - initialize TDT device for sequenced play
%   Y = SeqplayInit(Fsam) initializes the RX6 for sequenced play by
%   loading the circuit. Fsam is the (approx) sample rate in kHz.
%   SeqplayInit returns a struct containing the exact sample frequency
%   in kHz and other details of the sequenced play setup. 
%   
%   SeqplayInit(Fsam, Dev) uses device Dev. Default Dev is sys3defaultdev.
%   
%   SeqplayInit(Fsam, Dev, Postfix) loads rcx circuit [Dev 'seqPlay' Postfix]
%   instead of the default circuits 'RP2seqPlay' or 'RX6seqPlay'.
%   For instance,
%        SeqplayInit(Fsam, 'RX6', '-triggering') 
%   loads the 'rx6seqplay-triggering.rcx' circuit, which allows triggering
%   of spike times over digital inputs.
%   SeqplayInit(Fsam, Dev, '') is the same as SeqplayInit(Fsam, Dev).
%   Special values for Postfix are those starting with a '#'. They refer to
%   settings stored in the Seqplay setup file of the current computer. 
%            Postfix = '#Foo' 
%   is equivalent to 
%            Postfix = fromsetupfile('seqplay','fooCircuitPostfix')
%   Analogously,
%            Postfix = '#Foo/mydefault' 
%   is equivalent to 
%            Postfix = fromsetupfile('seqplay','foo', 'mydefault')
%   which applies the default 'mydefault' if the local setup value 'Foo' is 
%   not found on this computer.
%   Finally, '#noload' assumes that the right circuit has already been
%   loaded and refrains from loading a new one.
%
%   [Y,Recycled]=SeqplayInit(Fsam, Dev, PostFix, 1) does not re-load the 
%   circuit if it has already been loaded using the right sample rate, but 
%   re-uses the present circuit. Note that the sample rate must match exactly. 
%   The logical output argument Recycled indicates whether a loaded circuit
%   was indeed re-used. This information is important for timing calibration.
%   (See sys3loadCircuit, calibrateRX6).
%   
%   Type 'help Seqplay' to get more detailed instructions on sequenced playback.
%
%   See also SeqplayList, SeqplayGO, SeqplayHalt, SeqplayStatus, RX6Calibrate.



% SPinfo: persistent structure containing Seqplay information

% default device
if nargin<2, Dev = ''; end 
if nargin<3, Postfix = ''; end 
if nargin<4, Recycle = 0; end 

if isempty(Dev), [dum, Dev] = sys3dev(''); end 

if isequal('#noload', Postfix),
    % reverse engineer postfix from loaded circuit
    Postfix = sys3CircuitInfo(Dev, 'CircuitFile');
    [dum, Postfix] = fileparts(Postfix);
    [dum, Postfix] = strtok(Postfix, '-');
elseif ~isempty(Postfix) && isequal('#', Postfix(1)),
    Postfix = Words2cell(Postfix(2:end), '/');
    if isscalar(Postfix), DefValue = '';
    else, DefValue = Postfix{2};
    end
    Postfix = FromSetupFile('Seqplay', [Postfix{1} 'CircuitPostfix'], '-default', DefValue);
end

% reset seqplay info shared by seqplayXXX functions
private_seqPlayInfo('clear');

% load circuit
% different device types need different circuits
DoesStamping = 0;
switch upper(Dev(1:3)),
    case 'RX6',
        Circuitname = ['RX6seqPlay' Postfix];
        DoesStamping = 1; % RX6 can do time stamping
    case 'RP2',
        Circuitname = ['RP2seqPlay' Postfix];
    otherwise,
        error(['No SeqPlay circuit available for device ''' Dev '''.']);
end
%dbstack, Circuitname
Fsam_requested = Fsam;
CircuitDir = fullfile(fileparts(which(mfilename)),'RPvdS'); % .\RPvdS as seen from this mfile
FullCircuitName = fullfile(CircuitDir,[Circuitname '.rcx']);
[Fsam, CycUse, RCLD] = sys3loadCircuit(FullCircuitName, Dev, RX6sampleRate(Fsam), Recycle); % returns actual sample rate
DCoffset = FromSetupFile('Seqplay','RX6_anti_DC','-default', [0 0]); % DC offset to compensate offset of D/A converter output

% store info in persistent SPinfo
Status = 'initialized';
ID = SetRandState; % unique ID for DACprogress etc
private_seqPlayInfo(CollectInStruct(Dev, Circuitname, CircuitDir, Fsam_requested, Fsam, CycUse, DoesStamping, DCoffset, '-', Status, ID));

% Initialize circuit in such a way it won't run away
sys3halt(Dev);
sys3setpar(-1, 'EndSample', Dev);
sys3write([0 1 -1], 'SwitchListL', Dev, 0, 'I32'); 
sys3write([0 1 -1], 'SwitchListR', Dev, 0, 'I32'); 
sys3write(DCoffset(1)+zeros(1,1000), 'WaveformsL', Dev);
sys3write(DCoffset(2)+zeros(1,1000), 'WaveformsR', Dev);
sys3write(zeros(1,1000), 'OffsetListL', Dev);
sys3write(zeros(1,1000), 'OffsetListR', Dev);

% if DoesStamping, % time-stamp-related buffers
%     sys3write('zero', 'SubStamp', Dev);
%     sys3write('zero', 'SamStamp', Dev);
%     sys3setpar(0, 'CalibRelay', Dev);
% end

y = private_seqPlayInfo;





