function [Y, Tsnip, Snips] = cutFromWaveform(dx, Y, Xcut, Xwin, makeNan);
% cutFromWaveform - cut out selected parts from waveform
%    Y=cutFromWaveform(dx, Y, Xcut, [-Xpre Xpost], makeNan) removes from 
%    waveform Y the chunks centered at the X-values in array Xcut, and 
%    replaces them by the linear interpolation of their end points or, if 
%    makeNan is true, by Nans. Default makeNan is false. The X-values start 
%    at 0 and have steps dx.
%    The Kth cuts starts at Xcut(K)-XPre and ends at Xcut(K)+Xpost. 
%    If multiple cut intervals overlap, they are
%    merged into one interval. Cut intervals exceeding the edges are
%    replaced by "horizontal" patches having the value of the end that is
%    within the X range. Intervals that fall completely outside the X range
%    are ignored.
%
%    cutFromWaveform([dx x0], ...) lets the X-axis start at x0 (see
%    timeaxis).
%
%    [Y, Xsnip, Snips] = Y = cutFromWaveform(...) also returns the chunks 
%    removed from Y as columns of matrix Snips. Xsnip is the corresponding
%    X-axis running from -Xpre to Xpost. Values falling outside the range
%    are replaced by NaNs.
%    
%    See also timeaxis, INTERP1.

[makeNan] = arginDefaults('makeNan', 0);

needSnips = (nargout>2); 
% handle default values of dx and determine x0
if isequal([], dx),
    dx = [1 0];
elseif isscalar(dx),
    dx = [dx 0];
end
x0 = dx(2);
dx = dx(1);
Nsam = numel(Y);
SzY = size(Y);
Y = Y(:);

% convert X values to sample indices
x2i = @(x)round((x-x0)/dx);
[Xcut, isort] = sort(Xcut(:));
i0 = x2i(Xcut+Xwin(1)); % start of cuts
i1 = x2i(Xcut+Xwin(2)); % ends
Tsnip = linspace(Xwin(1), Xwin(2), round(diff(Xwin)/dx)+1).';
% save snips if requested
if needSnips, 
    Snips = [];
    for icut = 1:numel(i0),
        [ii0, ii1] = deal(i0(icut), i1(icut)); % formal start & end samples (may be out of range)
        iii0 = max(1,ii0); iii1 = min(ii1,Nsam); % true start & end samples (forced to fit)
        sn = [nan(iii0-ii0,1); Y(iii0:iii1); nan(ii1-iii1,1)];  % target chunks, completed w nans if out of range
        %dsize(Tsnip, Snips, sn)
        Snips = [Snips, sn];
    end
    % undo sorting
    sorti(isort) = 1:numel(isort);
    Snips = Snips(:,sorti);
end
ievent = VectorZip(i0,i1); % alternate starts and ends
while 1, % overlapping cut intervals - fuse them
    jbackward = find(diff(ievent)<=0);
    if isempty(jbackward), break; end
    ievent([jbackward(:); 1+jbackward(:)]) = [];
end
i0 = ievent(1:2:end);
i1 = ievent(2:2:end);
for icut = 1:numel(i0),
    [ii0, ii1] = deal(i0(icut), i1(icut));
    if ii0<1 && ii1>Nsam,
        error('Length of cut is larger than waveform from which to cut.');
    elseif ii0<1, % extrapolate end value
        if makeNan, Y(1:ii1) = nan;
        else, Y(1:ii1) = Y(ii1);
        end
    elseif ii1>Nsam, % extrapolate start value
        if makeNan, Y(ii0:end) = nan;
        else, Y(ii0:end) = Y(ii0);
        end
    else, % interpolate between end values
        if makeNan, Y(ii0:ii1) = nan;
        else, Y(ii0:ii1) = interp1([ii0,ii1], Y([ii0,ii1]), ii0:ii1);
        end
    end
end

Y = reshape(Y, SzY);










