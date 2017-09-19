function figh = BasicChunkPlot(Chunk, dt, sam, Ynoise, varargin);
% BasicChunkPlot - plot samples and highlight chunks
%    BasicChunkPlot(C, dt, S, Ynoise) plots the samples in S against time (ms).
%    The time axis starts at zero and has a spacing of dt ms. Horizontal
%    dashed lines, ranging form time 0 to the end, are plotted at vertical
%    positions Ynoise(k). The chunks are highlighted in red. Additional
%    plot args will be passed to the call to xplot that renders the chunks.
%    
%    BasicChunkPlot opens a new figure and returns its graphics handle.
%
%    ChunkPlot(C, [dt t0], S, Ynoise) uses a starting time t0 instead of
%    zero (see timeAxis).
%
%    Note that cutting the chunks from S (and, possibly, T) is based on the
%    sample index fields of C, C.istart and C.iend. Any time information
%    contained in C is ignored.
%
%    See also ChunkPlot, xplot, TimeAxis, Xchunks, GetChunk.

figh = figure;
set(figh,'units', 'normalized', 'position', [0.0328 0.49 0.948 0.381]);
sam = sam(:);
T = timeaxis(sam,dt);
plot(T, sam);
% dashed horizontal lines
for iline=1:numel(Ynoise),
    xplot(xlim, [0 0]+Ynoise(iline), 'k:');
end
% Collect chunk waveforms and glue them together. Use nans as fast
% separators
X = []; Y = [];
for ic=1:numel(Chunk),
    i0 = Chunk(ic).istart;
    i1 = Chunk(ic).iend;
    X = [X; T(i0:i1); nan];
    Y = [Y; sam(i0:i1); nan];
end
xplot(X,Y,'r', varargin{:}); ylim(ylim);
xlabel('Time (ms)');
