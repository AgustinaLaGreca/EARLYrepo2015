function seqplayPlot(flag);
% SeqplayPlot - plot buffer contents of seq play (debug function)
%
% SeqplayPlot('played') attempts to retrieve the samples visited (RX6 only)


if nargin<1, flag=''; end
SPinfo = seqplayinfo;

dt = 1/SPinfo.Fsam; % sample period in ms

played = isequal('played', flag);
if played, error NYI, end

% Left
nsamL = SPinfo.NsamL;
if ~isempty(nsamL),
    offsetL = SPinfo.WavL_offset;
    % read L waveform
    WAV = []; NSAM = [];
    for iplay=1:length(SPinfo.iwaveL),
        iwave = SPinfo.iwaveL(iplay);
        nrep = SPinfo.NrepL(iplay);
        offset = SPinfo.WavL_offset(iwave);
        nsam = SPinfo.NsamL(iwave);
        NSAM = [NSAM repmat(nsam,1,nrep)];
        wav = sys3read('WaveformsL', nsam, SPinfo.Dev, offset, 'F32');
        WAV = [WAV, repmat(wav,1,nrep)];
    end
    istart = [0 cumsum(NSAM(1:end-1))]; % start indices of waveform buffers
    dplot(dt,WAV);
    xplot(dt*istart, WAV(istart+1), 'o');
end

% Rite
nsamR = SPinfo.NsamR;
if ~isempty(nsamR),
    offsetR = SPinfo.WavR_offset;
    % read L waveform
    WAV = []; NSAM = [];
    for iplay=1:length(SPinfo.iwaveR),
        iwave = SPinfo.iwaveR(iplay);
        nrep = SPinfo.NrepR(iplay);
        offset = SPinfo.WavR_offset(iwave);
        nsam = SPinfo.NsamR(iwave);
        NSAM = [NSAM repmat(nsam,1,nrep)];
        wav = sys3read('WaveformsR', nsam, SPinfo.Dev, offset, 'F32');
        WAV = [WAV, repmat(wav,1,nrep)];
    end
    istart = [0 cumsum(NSAM(1:end-1))]; % start indices of waveform buffers
    dplot(dt,WAV,'r');
    xplot(dt*istart, WAV(istart+1), 'or');
end


