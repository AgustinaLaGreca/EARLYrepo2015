function C = EvalChunks(C, S, FN, fhan, varargin);
% EvalChunks - evaluate chunks of a signal with respect to given property
%    C = EvalChunks(C, S, 'Foo', @FCN, ...) evaluates chunks of array S as
%    specified by struct arrray C. The k-th chunk of S consists of the 
%    segment S(C(k).istart:C(k).iend). A new field Foo is added to C, whose
%    values are a chunk-wise evaluation of function FCN. The first argument
%    passed to FCN is the segment, the second is the single chunk C(k).
%    Optional arguments indicated by ... are passed also to FCN.
%    Special values are the strings 'PrevChunk' and 'NextChunk' for which
%    C(k-1) and C(k+1) are substituted, respectively, or [] for k=1 and
%    k=numel(C), respectively. Note that the chunks need not be adjacent.
%
%    C = EvalChunks(C, S, {'Fd1' 'Fd2' .. 'FdN'}, @FCN, ...) creates
%    multiple new fields Fd1, Fd1.. FdN by calling function FCN, which
%    should provide N output arguments.
%
%    See also Xchunks, ChunkPlot, Getchunk.

if ~isfhandle(fhan),
    error('Fourth input arg must be function handle.');
end

% look for special values in varargin (see help text)
iprev = [] ; inext = [];
for ii=1:numel(varargin),
    if isequal('PrevChunk', varargin{ii}), iprev=[iprev, ii];
    elseif isequal('NextChunk', varargin{ii}), inext=[inext, ii]; end
end
anyPrev = ~isempty(iprev);
anyNext = ~isempty(inext);

[FN, Nf] = cellify(FN);

Nch = numel(C);
FV = cell(Nf,Nch); % the new field values end up in FV
for ic=1:Nch,
    seg = getchunk(S,C(ic));
    if anyPrev,
        if ic==1, Prv = []; else, Prv = C(ic-1); end
        [varargin{iprev}] = deal(Prv);
    end
    if anyNext,
        if ic==Nch, Nxt = []; else, Nxt = C(ic+1); end
        [varargin{inext}] = deal(Nxt);
    end
    [FV{:, ic}] = fhan(seg, C(ic), varargin{:});
end
% put values of FV in FN field of C
for ifid=1:Nf,
    [C.(FN{ifid})] = deal(FV{ifid, :});
end


