function [Chunk, figh] = xchunks(sam, Nsmooth, M, dt);
% xchunks - basic chunk extraction from recording trace
%   CH = xchunks(D, Nsmooth, M) extracts "chunk" array CH from data in array X. 
%   A chunk is a segment of the recording that is monotonically rising or 
%   falling. Not all of such monotonic segments end up as chunks; in order
%   to become a chunk, a segment must meet one of the following criteria:
%
%     1. its "span" (absolute value between start and end value) exceeds M
%        times the IQR of the recording;
%     2. The span of both "neighbors" (segments adjacent to it) exceeds M
%        times the IQR.
%     3. the segment is one of a pair of neighbors, both of which are "small"
%        (fail to meet criterium #1), but whos summed span exceed M*IQR.
%
%   Input args:
%          X: array of samples (typically a recording).
%    Nsmooth: length of smooting hamming window. 
%          M: IQR factor as described above. Default M=5. If M is a
%             2-component array, M(1) and M(2) are used as factors for
%             testing criterium 1&2 and 3, respectively.
%    Output arg:
%        CH: struct array of chunks.
%  
%    Xchunk performs the generic feature extraction that provides a basis 
%    of more specific feature extraction in terms of physiological
%    interpretation.
%
%    [CH, h] = Xchunks(X, Nsmooth, M, dt) also plots the waveform X using 
%    a sample period of dt ms. The chunks are highlited in red. The
%    graphics handle h of the newly opened figure is returned. The plotting
%    is delegated to BasicChunkPlot.
%
%    See also Smoothen, EvalChunks, GetChunk, BasicChunkPlot, ChunkPlot.

if nargin<3, M = 5; end % default IQR factor
if nargin<4, dt = []; figh=[]; end % default: no plotting

sam = sam(:);

MeanSam = mean(sam);
MedSam = median(sam); % median value of samples, serving as "baseline"
IqrSam = iqr(sam); % inter quartile range, serving as "natural spread"
Nsam = numel(sam);

% the threshold for the two criteria described in the help text
M
THR1 = M(1)*IqrSam
THR2 = M(end)*IqrSam

% chop up wave into rising and falling chunks. 
smsam = smoothen(sam, Nsmooth); % smoothing is needed to avoid senseless small chunks.
[imin ymin] = localmax(1,-smsam);
[imax ymax] = localmax(1,smsam);
clear smsam; % save some memory
% sort all extrema according to their time order
[yturn, iturn] = deal([-ymin(:).', ymax(:).'], [imin(:).', imax(:).']);%all extrema + sample indices in row vector
[yturn, iturn] = sortAccord(yturn, iturn, iturn); % sort according to time
% Candidate chunks are the segments in between the turning points. 
%f1; plot(sam); xplot(iturn, yturn, 'r*'); 
span = abs(diff(yturn)); % the span of each of these segments.
pass1 = (span>THR1); % segments having big span meet first criterion
pass2 = [false, (pass1(1:end-2) & pass1(3:end)), false]; % having two big neighbors makes you a chunk, too
span(pass1) = 0; % the effects of big segments are exhausted. Discard their bookkeeping.
dbspan = span(1:end-1)+span(2:end); % "double span": sum over pairs of adjacent segments; big segments excluded
%f2; hist(span/IqrSam, 100);
% test 2nd criterium: adjacent small segments whose joint span exceeds THR2 both count as chunks.    
pass3 = [(dbspan>THR2), false] | [false, (dbspan>THR2)]; % either segment contributing to a large joint span 
% Initialize chunk array
ichunk = find(pass1 | pass2 | pass3); % indices of chunks
%ichunk = find(pass1 | pass2);
Chunk = struct('istart', num2cell(iturn(ichunk)), ...
    'iend', num2cell(iturn(ichunk+1))  );

if ~isempty(dt),
    BaseLines = MedSam+IqrSam*[-1 1]*M(1)/2;
    figh = BasicChunkPlot(Chunk, dt, sam, BaseLines);
end








