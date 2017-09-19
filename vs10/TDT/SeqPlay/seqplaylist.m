function SPinfo = SeqplayList(iwaveL, NrepL, iwaveR, NrepR);
% SeqplayList - specify play list for sequenced play
%   SeqplayPlaylist(iwaveL, nrepL, iwaveR, NrepR) specifies a play list with:
%     iwavL:  vector containing indices>=1 indicating the left-channel
%             waveforms previously loaded by SeqplayUpload. A single
%             index may occur multiple times in the list.
%     nrepL:  respective repetition counts for iwavL.
%     iwavR:  index vector for right channel.
%     nrepR:  repetition vector for right channel.
% 
%   To play over the right channel only, leave both iwavL and nrepL empty.
%   To play over the left channel only, leave both iwavR and nrepR empty,
%   or omit these input arguments.
%
%   For a play list to be specified, the buffers must have been uploaded. 
%   The same collection of buffers can be used with different play lists 
%   without having to upload the waveforms again; a call to SeqplayList 
%   overrides any previous list specification.
% 
%   Type 'help Seqplay' to get an overview of sequenced playback.
%   Sing is an audio test that uses SeqPlay functionality.
%
%   Example
%
%      SeqPlayList([],[],[1 2 1], [2 10 1]); 
%
%   This call activates only the right channel. Waveform #1 is first played
%   twice. Next Waveform #2 is played 10 times. Finally Waveform #1 is
%   played once.
% 
%   See also Seqplay, Sing, SeqPlayInit, SeqplayStatus.

if nargin <3, 
    iwaveR = [];
    NrepR = [];
end

SPinfo = private_seqPlayInfo; % info shared by seqplayXXX
if isempty(SPinfo), 
    error('Sequenced play not initialized. Call SeqPlayInit first.');
elseif isequal('initialized', SPinfo.Status),
    error('No waveforms loaded. Call SeqplayUpload first.');
end

seqplayhalt; % abort any ongoing D/A and make sure the circuit will not run away due to alterations of play lists

if isempty(iwaveL) && isempty(iwaveR),
    error('Play lists for both channels are empty.');
end
% real work is delegated to local function
[NsamTotPlayR, iwaveR_full, itickSwitchR, iplayR, irepR] = local_doit(iwaveR, NrepR, 'R', SPinfo);
[NsamTotPlayL, iwaveL_full, itickSwitchL, iplayL, irepL] = local_doit(iwaveL, NrepL, 'L', SPinfo);

% if [] iwas specified for one channnel, it should output a sufficient ...
% ... number of zeros - formally that is waveform #0
if isempty(iwaveL),
    Nrep = ceil(NsamTotPlayR/SPinfo.Nzero); % number of times the zero buffer shoud be repeated to match duration of play
    [NsamTotPlayL, iwaveL_full, itickSwitchL] = local_doit(0, Nrep, 'L', SPinfo);
    NsamTotPlayL = NsamTotPlayR; % don't worry - if anything there are too *many* zeros in L
elseif isempty(iwaveR), % same thing with channels swapped
    Nrep = ceil(NsamTotPlayL/SPinfo.Nzero);
    [NsamTotPlayR, iwaveR_full, itickSwitchR] = local_doit(0, Nrep, 'R', SPinfo);
    NsamTotPlayR = NsamTotPlayL; % don't worry - if anything there are too *many* zeros in L
end

if ~isequal(NsamTotPlayL, NsamTotPlayR),
    error('Unequal total # samples in left & right channels.');
end

c____________ = '_______________';
d____________ = '_______________';
NsamTotPlay = NsamTotPlayL;
private_seqPlayInfo(CollectInStruct(NsamTotPlay, ...
    c____________, ...
    iwaveL, NrepL, iwaveL_full, iplayL, irepL, itickSwitchL, ...
    d____________, ...
    iwaveR, NrepR, iwaveR_full, iplayR, irepR, itickSwitchR));
private_seqPlayInfo('Status', 'listed');

sys3run(SPinfo.Dev);

SPinfo = private_seqPlayInfo;
%----------------------------
function [NsamTotPlay, iwave, itick, iplay, irep] = local_doit(iwave, Nrep, Chan, SPinfo);
if isempty(iwave), [NsamTotPlay, iwave, itick, iplay, irep] = deal([]); return; end; % postpone silent D/A unitil other channel is known
% assemble seq play instruction & upload to circuit
if length(iwave)~=length(Nrep),
    error('Wave index vector and Nrep vector have unequal length.');
end

% expand the reps into single buffer playbacks
[iwave, iplay, irep] = local_expand_reps(iwave, Nrep);
Nsam = SPinfo.(['Nsam' Chan]); % NsamL or NsamR
Nsam = [SPinfo.Nzero Nsam]; % treat heading silence as "zeroth buffer"
if any(iwave>length(Nsam)), error('Elements of iwave exceed # stored waveforms'); end;
% # samples of these single plays
Nsam = Nsam(1+iwave); % 1+ corrects for prepending to Nsam of SPinfo.Nzero above
% offset of the respective buffers these single plays
offsets = SPinfo.(['Wav' Chan '_offset']); % WavL_offset or WavR_offset 
offsets = [0 offsets];  % as before: treat heading silence as "zeroth buffer"
offsets = [offsets(1+iwave) 0]; % 1+ corrects for prepending in previous line; the final zero is the return to silence @ the end
% compute the tick values at the start of each of the single play (at the "jumps")
itick = cumsum([2 Nsam]);
% the Master Clock of the circuit is not reset on jump events; the shift ...
% ... from tick to isample must take this into account
ishift = offsets-itick;
% tell the circuit
iswitch = [0 -1 itick-1]; % first 2 values handle spurious updating @ startup
ishift = [0 0 0 ishift];
sys3write(iswitch, ['SwitchList' Chan], SPinfo.Dev, 0, 'I32');
sys3write(ishift, ['OffsetList' Chan], SPinfo.Dev, 0, 'I32');
if isequal('RX6seqPlay-divided', SPinfo.Circuitname), % HACK: multi-processor data transfer costs time; provide that time to avoid click artefacts
    for ii=1:round(600/SPinfo.Fsam),
        iswitchCheck = sys3read(['SwitchList' Chan], length(iswitch), SPinfo.Dev, 0, 'I32');
    end
end
NsamTotPlay = 2+sum(Nsam);

function [IWAVE, iplay, irep] = local_expand_reps(iwave, nrep);
% expand reps using ones (NOTE v*ones(..) is faster than repmat(v,..))
qq = arrayfun(@(iw,nr)(iw*ones(1,nr)), iwave, nrep, 'UniformOutput', false); IWAVE = [qq{:}];
qq = arrayfun(@(iw,nr)(iw*ones(1,nr)), 1:numel(iwave), nrep, 'UniformOutput', false); iplay = [qq{:}];
qq = arrayfun(@(nr)(1:nr), nrep, 'UniformOutput', false); irep = [qq{:}];













