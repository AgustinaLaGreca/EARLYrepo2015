function SPinfo = SeqplayUpload(LBuff, RBuff, scaleL, scaleR);
% SeqplayUpload - Upload samples for sequenced play over TDT device
%   SeqplayUpload({SL1,SL2,...}, {SR1,SR2,...}) uploads waveforms SL1, SL2, .. 
%   to the left channel, and waveforms SR1, SR2, .. to the right channel, in
%   preparation for sequenced playback. Individual waveforms must be column vectors
%   containing at least 2 samples. Different waveforms may have different lengths.
%   The samples are uploaded to the device.
%   SeqplayUpload returns a struct containing details of the sequenced play setup.
%
%   SeqplayUpload(SL, SR, scaleL, scaleR) first scales the samples before 
%   uploading them. Waveform SL{k} is multiplied by scaleL(k) and waveform
%   SR{k} is multiplied by scaleR(k).
%
%   For monaural playback, upload an empty buffer cell {} to the non-active channel.
%
%   See also SeqPlay, SeqplayInit, Sing, SeqplayList, SeqplayGO, SeqplayHalt, SeqplayStatus.

if nargin<3, scaleL=1; end
if nargin<4, scaleR=1; end

if ~isequal( 'cell', class(LBuff) ) | ~isequal( 'cell', class(LBuff) ),
   error('Input arguments should be cell arrays.');
end

Nzero = 1000; % # samples of silence buffer

SPinfo = private_seqPlayInfo; % info shared by seqplayXXX
if isempty(SPinfo), error('Sequenced play not initialized. Call SeqPlayInit first.'); end
seqplayhalt; 

LeftEmpty = isempty(LBuff);
RightEmpty = isempty(RBuff);
if LeftEmpty & RightEmpty,
   error('Both sides are empty. Upload aborted.');
end
% sys3setpar(double(~LeftEmpty), 'ScaleL', SPinfo.Dev);
% sys3setpar(double(~RightEmpty), 'ScaleR', SPinfo.Dev);

% ---check waveforms 
if any(~cellfun(@isvector, LBuff)), error('Waveform buffers must be numerical column vectors.'); end
if any(~cellfun(@isvector, RBuff)), error('Waveform buffers must be numerical column vectors.'); end
NsamL = cellfun(@numel, LBuff); % vector containing respective length of left buffers
NsamR = cellfun(@numel, RBuff); % idem right
NwavL = numel(LBuff); NwavR = numel(RBuff);
if any(NsamL<2) || any(NsamR<2), error('Each waveform buffer must contain at least 2 samples.'); end
WavL = real(cat(1, zeros(Nzero,1), LBuff{:})); % grand left waveform; ignore any imag part
WavR = real(cat(1, zeros(Nzero,1), RBuff{:})); % grand right waveform; ignore any imag part
% if (length(WavL) > 1e6) || (length(WavR) > 1e6), 
%    error('Too many samples. Max = 1e+006-Nzero per channel.');
% end
% apply segment-wise scaling if requested
WavL = local_Scale(WavL, NsamL, scaleL, Nzero); 
WavR = local_Scale(WavR, NsamR, scaleR, Nzero); 
%dsize(WavL, WavR);
% ---store waveforms in circuit
sys3write(SPinfo.DCoffset(1)+WavL, 'WaveformsL', SPinfo.Dev, 0, 'F32', 3);  % 3= check for corruption; on failure ..
sys3write(SPinfo.DCoffset(2)+WavR, 'WaveformsR', SPinfo.Dev, 0, 'F32', 3); % .. ; try at most 3 times.
% silent channels: make sure no junk is left
% LeftEmpty, RightEmpty
if LeftEmpty, sys3write(SPinfo.DCoffset(1)+zeros(2*Nzero,1), 'WaveformsL', SPinfo.Dev, 0, 'F32', 3); end
if RightEmpty, sys3write(SPinfo.DCoffset(2)+zeros(2*Nzero,1), 'WaveformsR', SPinfo.Dev, 0, 'F32', 3); end


% ---store SeqPlay info & return it
WavL_offset = cumsum([Nzero NsamL]); % start indices
WavR_offset = cumsum([Nzero NsamR]); 
Status = 'loaded';
b____________ = '_______________';
private_seqPlayInfo(CollectInStruct(NsamL, NsamR, LeftEmpty, RightEmpty, Nzero, WavL_offset, WavR_offset, Status, b____________));
SPinfo = private_seqPlayInfo;

%======LOCALS======================================
function W = local_Scale(W, Nsam, Scale, Nzero); % apply segment-wise scaling if requested
Nwav = numel(Nsam);
if ~isequal(1, Scale),
    if Nwav==1,
        W = W*Scale;
    elseif Nwav>0, % different buffer segments have different scale factor
        scaleOffset = cumsum([Nzero Nsam(1:end-1)]);
        for ii=1:Nwav,
            % scaleOffset(ii)+[1 Nsam(ii)]
            ir = scaleOffset(ii)+(1:Nsam(ii)); % sample range
            W(ir) = W(ir)*Scale(ii);
        end
    end
end

 





