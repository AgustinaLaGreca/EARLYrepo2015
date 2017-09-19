function S = seqplayStatus(Chan, Prop);
% SeqplayStatus - current status of sequenced play
%   S = SeqplayStatus(Chan) returns information on the
%   current status of the SeqPlay in struct S.
%
%        Active: indicates if playback is active or not.
%     WaveCount: current waveform counted in order of play
%     WaveIndex: buffer index of current waveform (cf seqplayupload)
%     IplayList: current item of playlist
%          Irep: repetition count (cf seqplaylist)
%      isam_abs: 0-based sample count since start-of-play
%      isam_rel: sample count re start of current waveform
%
%   S = SeqplayStatus(Chan, 'Foo') anly returns field Foo.
%
%   See also SeqPlay, SeqPlayGo, SeqPlayDAwait.

[Chan, Prop] = arginDefaults('Chan/Prop', 1 ,'');
    
if isequal('L', Chan), Chan=1;
elseif isequal('R', Chan), Chan=2;
else, Chan = 1:2;
end

[Active] = 0; % default: not playing anything we know of
[WaveCount, WaveIndex, IplayList, Irep, isam_rel] = deal([0 0]); % defaults of binaural params
isam_abs = 0; % not this one..
S = CollectInStruct(Active, isam_abs, WaveCount, WaveIndex, IplayList, Irep, isam_rel);

SPinfo = private_seqPlayInfo; % info shared by seqplayXXX
if isempty(SPinfo), 
    if ~isempty(Prop), S = [S.(Prop)]; end
    return; 
end; % not initialized at all
% don't believe SPinfo, but check whether the circuit is loaded
CircuitInfo = sys3CircuitInfo(SPinfo.Dev);
if isempty(CircuitInfo), % circuit was cleared behind our back
    private_seqPlayInfo('clear'); 
    SPinfo=[]; 
    if ~isempty(Prop), S = [S.(Prop)]; end
    return; % circuit
elseif ~isequal('listed', SPinfo.Status), 
    if ~isempty(Prop), S = [S.(Prop)]; end
    return; % not ready for D/A
end

% potentially active: find out where we are timewise

isam_abs = sys3getpar('Tick', SPinfo.Dev); % absolute time in samples since start of D/A

Active = ~isequal(0, isam_abs) && (SPinfo.NsamTotPlay+1>isam_abs); % between end and begin sample
if ~Active, S = CollectInStruct(Active, WaveCount, WaveIndex, IplayList, Irep, isam_abs, isam_rel); 
    if ~isempty(Prop), S = [S.(Prop)]; end
    return; 
end

% active: find out where we are waveform wise
iL = find(SPinfo.itickSwitchL<=isam_abs, 1, 'last'); if isempty(iL), iL=1; end
iR = find(SPinfo.itickSwitchR<=isam_abs, 1, 'last'); if isempty(iR), iR=1; end
WaveCount = [iL iR];
Nwave = [length(SPinfo.iwaveL_full) length(SPinfo.iwaveR_full)];
WaveCount = min([iL iR], Nwave);
iL  = WaveCount(1); iR  = WaveCount(2);
% now identify what is being played, where that is in the playlist, and
% which rep it is.
WaveIndex = [SPinfo.iwaveL_full(WaveCount(1)) SPinfo.iwaveR_full(WaveCount(2))];
if isempty(SPinfo.iwaveL), % silent L channel
    IplayListL = 0;
    IrepL = 0;
    isam_relL = 0;
else,
    IplayListL = SPinfo.iplayL(iL);
    IrepL = SPinfo.irepL(iL);
    iwa = WaveIndex(1);
    if iwa==0, % zeroth buffer has length Nzero
        Nsam = SPinfo.Nzero;
    else,  % # samples in currrent waveform
        Nsam = SPinfo.NsamL(iwa); 
    end
    isam_relL = rem(isam_abs-SPinfo.itickSwitchL(iL), Nsam);
end
if isempty(SPinfo.iwaveR), % silent R channel
    IplayListR = 0;
    IrepR = 0;
    isam_relR = 0;
else,
    IplayListR = SPinfo.iplayR(iR);
    IrepR = SPinfo.irepR(iR);
    iwa = WaveIndex(2);
    if iwa==0, % zeroth buffer has length Nzero
        Nsam = SPinfo.Nzero;
    else,  % # samples in currrent waveform
        Nsam = SPinfo.NsamR(iwa); 
    end
    isam_relR = rem(isam_abs-SPinfo.itickSwitchR(iR), Nsam);
end
IplayList = [IplayListL, IplayListR]; 
Irep = [IrepL, IrepR];
isam_rel = [isam_relL, isam_relR];
WaveCount = WaveCount(Chan); WaveIndex = WaveIndex(Chan);
IplayList = IplayList(Chan); Irep = Irep(Chan);
isam_rel = isam_rel(Chan);

NsamTot=SPinfo.NsamTotPlay;
dt = 1/SPinfo.Fsam; % sample period in ms
S = CollectInStruct(Active, WaveCount, WaveIndex, IplayList, Irep, isam_abs, isam_rel, Chan, NsamTot, dt);

if ~isempty(Prop), S = [S.(Prop)]; end




