function [Tpeak, CV] = findtypsnip(dt, Y, Width);
% scramblePhase - scrambl phase of time signal.
%    T = findtypsnip(dt,Y, Width) attempts to find a typical snippet T
%    for describing waveform Y. 
%    dt is sample period of Y; Width is duration of T.
%
%    See scramblePhase.

Nit = 500;

Nsam = numel(Y);
NsamT = round(Width/dt);
Ys = scramblePhase(Y);

istart = RandomInt(Nsam-NsamT, [1 Nit]);
for ii=1:Nit,
    Ti = local_Ti(Ys, istart(ii), NsamT);
    % use Ti as filter for Y and evaluate the peaks
    [Npeak(ii), MeanPeakSize(ii), MeanBigPeakSize(ii)] = local_CV(dt, Y, Ti, Width);
end

set(figure,'units', 'normalized', 'position', [0.606 0.146 0.402 0.468])
lh = plot(Npeak, MeanBigPeakSize, '.');
IDpoints(lh, {Ys NsamT dt}, istart, @(Ys,N,dt,i)['#' num2str(ii)], 'show snippet', @local_plot)

% find "best" snippet & plot it
[dum imaxPz] = max(MeanBigPeakSize)
local_plot(Ys, NsamT, dt, istart(imaxPz), 'k', 'linewidth', 2);
[Npeak, MeanPeakSize, MeanBigPeakSize, Tpeak, CV] = local_CV(dt, Y, Ti, Width);
[dum, imax] = max(Ti);
Tpeak = Tpeak-dt*imax;

% ==================
function Ti = local_Ti(Ys, istart, NsamT);
iend = istart+NsamT-1;
Ti = Ys(istart:iend); % candidate typical snip
% normalize Ti to make competition fair
Ti = Ti-mean(Ti);
Ti = Ti/std(Ti);

function [Npeak, MeanPeakSize, MeanBigPeakSize, Tpeak, CV] = local_CV(dt, Y, Ti, Width);
CV = conv(Y,Ti);
CV = smoothen(CV,Width/7, dt);
%dsize(CV,dt,Y, Ti)
[Tpeak, Ypeak]=localmax(dt, CV);
Npeak = numel(Tpeak);
MeanPeakSize = mean(Ypeak);
%hist(Ypeak,50); pause; clf
SortedPeaks = sort(Ypeak, 'descend');
MeanBigPeakSize = mean(SortedPeaks(1:round(Npeak/4)).^2);

function local_plot(Ys, NsamT, dt, istart, varargin);
%dsize(Ys, NsamT, dt, istart)
hf = findobj('tag', 'SnippetS');
if isempty(hf), 
    set(figure, 'tag', 'SnippetS'); 
    set(gcf,'units', 'normalized', 'position', [0.229 0.213 0.359 0.322]);
else, 
    figure(hf); 
end
Ti = local_Ti(Ys, istart, NsamT);
if isempty(varargin),
    xdplot(dt, Ti(end:-1:1), lico(istart));
else,
    xdplot(dt, Ti(end:-1:1), varargin{:});
end










