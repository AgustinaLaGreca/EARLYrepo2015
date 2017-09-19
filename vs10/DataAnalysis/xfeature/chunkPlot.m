function C = chunkPlot(C, T, S, varargin);
% ChunkPlot - plot or highlight chunks
%    ChunkPlot(C, T, S, ...) xplots chunks in a S-versus-T plot. 
%    S is the Signal from which the chunks are cut using struct array C 
%    (see Xchunks, Getchunk). T is typically the time axis of recorded 
%    signal S. If T is a single number or a 2-element vector, a time axis
%    TimeAxis(S,T) is used. Alternatively T and S may both be signals, 
%    e.g., a recording and its time derivative. When S and/or T are 
%    matrices, the chunk extraction is done along the first dimension.
%
%    Note that cutting the chunks from S (and, possibly, T) is based on the
%    sample index fields of C, C.istart and C.iend. Any time information 
%    contained in C is ignored. 
%
%    Trailing args indicated by ellipses (...) are passed to xplot.
%
%    See also xplot, Xchunks, GetChunk, EvalChunks.

if numel(T)<=2 && numel(T)<numel(S), % single-number T means sample period; pair means [dt t0] (see timeaxis)
    T = timeaxis(S, T); % expand to full time axis
end
% make sure that time dimension is 1st dim
if isvector(T), T=T(:); end
if isvector(S), S=S(:); end
% use Nans as "visual "separators" to avoid calling PLOT for each chunk.
TNaNa = nan(1,size(T,2));
SNaNa = nan(1,size(S,2));

Nch = numel(C);
TT = []; SS = [];
for ic=1:Nch,
    c = C(ic);
    Tseg = getchunk(T, c);
    Sseg = getchunk(S, c);
    TT = [TT; Tseg; TNaNa];
    SS = [SS; Sseg; SNaNa];
end
%dsize(TT,SS), varargin{:}
xplot(TT, SS, varargin{:});





