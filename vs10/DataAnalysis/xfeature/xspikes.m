function S = xspikes(D, isweep, ichan, Vmax, doPlot, AlignThr);
% Xspikes - extract spikes from recorded traces
%   xspikes(D, isweep, ichan), where D is return argument of readABF, 
%   attemps to extract complex spikes from the specified sweep & channel.
%   Default isweep=1; ichan=1.
%
%   xspikes(D, isweep, ichan, Vmax, doPlot, AlignThr) restricts the 
%   analysis to events whose size exceeds Vmax. Default Vmax is 
%   0.1*max(abs(hilbert(Y))), where Y is the raw time signal.
%   doPlot specifies whether the results are plotted.
%   AlignThr is threshold criterion for time-alignment. Default is
%   0.5*min(Vmax).
%
%   See also readABFdata, ABFplot.

if nargin<2, isweep=1; end
if nargin<3, ichan=1; end
if nargin<4, Vmax=[]; end
if nargin<5, doPlot=1; end % default: do plot
if nargin<6, AlignThr=[]; end;

if numel(findobj('type','figure'))>15, warning('deleting figures!'); aa; end
% collect set of params that makes the calculation unique. Need to
% distinguish between data read by readTKABF versus readABFdata; they may
% stem from the same ABFfile, but the former is timealigned and scaled; the
% latter is not.
CacheParam = {D.Header.XX_FullFileName, D.ExpID, D.icond, D.NsamPerSweep, isweep, ichan, Vmax, AlignThr}; 
CFN = mfilename; % name of cache file
S = getcache(CFN, CacheParam);
if ~isempty(S), 
    if doPlot, xspikeplot(S); end;
    return; 
end

Tramp = 0.2; % ms ramps for Hilbert transform
maxLag = 1.5; % ms range for autocorr
Wfactor = 3; % width factor for smoothing of abs(y)
CorCrit = 0.9; % corr criterion for testing events

ad = [];
for iD=1:numel(D),
    ad = [ad; D(iD).AD(isweep, ichan)]; % the requested sweep(s)
end
D = D(1);
y = cat(1,ad.samples); % recorded potential in mV
dt = D.dt_ms; % time spacing in ms
Nsam = numel(y);
t = timeaxis(y, dt, 0);

% subtract offset, provide ramps and take Hilbert transform
Ymean = mean(y);
y = y-Ymean;
Nramp = 1+round(Tramp/dt);
tr = sin(2*pi*(1:Nramp).'/Nramp/4).^2;
%dplot(dt,tr, '.');
y(1:Nramp) = tr.*y(1:Nramp);
y(end+1-(1:Nramp)) = tr(Nramp:-1:1).*y(end+1-(1:Nramp));
y = hilbert(y);
yabs = abs(y);

if isempty(Vmax), Vmax = [0.1 0.5]*max(yabs); end % see help

% plot traces with "envelope" and Hilbert trf
%figure; hf1=gcf; set(gcf,'units', 'normalized', 'position', [0.026563 0.32422 0.51797 0.5752]);
% xplot(t, real(y), 'linewidth', 2);
% xplot(t, imag(y),'g:');
% xplot(t, yabs,'r');
% xplot(t, -yabs,'r');

% get complex autocorr fnc and plot its abs value; estimate typical spike
% width from it
figure; hf2=gcf; set(gcf,'units', 'normalized', 'position', [0.81094 0.13086 0.22813 0.19434]);
imaxlag = round(maxLag/dt);
[AC, ilag] = xcorr(y, imaxlag, 'coeff');
tlag = dt*ilag;
xplot(tlag, abs(AC));
W = abs(AC); % weights for averaging the time
SpikeWidth = sqrt(sum(tlag.^2*W)/sum(W));

% use spike width to smooth the absolute value of the recordings, and use
% this to time "events"
Wsm = exp(-0.5*(Wfactor*tlag/SpikeWidth).^2);  
xplot(tlag,Wsm,'r')
Wsm = Wsm/sum(Wsm); % normalize ..
ysmooth = conv(yabs,Wsm);
tsmooth = dt*(0:numel(ysmooth)-1).'-maxLag;
%figure(hf1); xplot(tsmooth, ysmooth, 'k');
[tmax, ymax, isort] = localmax(tsmooth, ysmooth);
%figure(hf1); xplot(tmax, ymax, 'sr')
figure; hf3=gcf; hist(ymax,100); set(gcf,'units', 'normalized', 'position', [0.51172 0.12793 0.30938 0.19336])
y_at_max = interp1(t,real(y),tmax);
%figure(hf1); xplot(tmax,y_at_max,'h')

% sort events to size; extract "snippets" of recording around events
tmax= tmax(isort); ymax= ymax(isort);
Nevent=numel(tmax);
%f4;set(gcf,'units', 'normalized', 'position', [0.39375 0.33398 0.58516 0.5791]);
iploco=1;
SnipRange = 2.5*SpikeWidth; % ms max time range around event center
NsamSnip = 1+2*round(SnipRange/dt);
t0snip = linspace(-SnipRange,SnipRange,NsamSnip);
Wsn = exp(-0.5*(t0snip/SpikeWidth).^2);
figure(hf2); xplot(t0snip, Wsn,'g'); Wsn = Wsn/sum(Wsn);
Npause = 10;
figure(hf3); YYY=[1 1]*mean(ylim);
icat=1;
Template{1} = [];
Ncat = numel(Vmax);

EventCount = 0; catIdx = zeros(Nevent,1); maxievent = 0;
for ievent = 1:Nevent
    if ymax(ievent)<Vmax(end+1-icat), 
        icat = icat+1;
        if icat>Ncat, maxievent = ievent-1; break; end
        Template{icat} = [];
        EventCount(icat) = 0;
    end
    col = ploco(icat);
    EventCount(icat) = EventCount(icat)+1;
    tsnip = tmax(ievent) + t0snip;
    ysnip = interp1(t,real(y),tsnip,'linear',0);
    %ysnip = ysnip-sum(ysnip.*Wsn);
    Template{icat} = [Template{icat}; ysnip];
    catIdx(ievent) = icat;
    %if ievent>1, f3; xplot(ymax(ievent+[-1 0]), YYY, col, 'linewidth', 4); end
    %f4; xplot(t0snip, ysnip, col);
end
tmax = tmax(1:maxievent);
ymax = ymax(1:maxievent);
catIdx = catIdx(1:maxievent); 
% compute median snippet for each magnitude class
% f5; set(gcf,'units', 'normalized', 'position', [0.0015625 0.3925 0.4375 0.525]);
for icat=1:Ncat,
    if EventCount(icat)>0,
        Template{icat} = median(Template{icat},1);
        %f5; xplot(t0snip+tshift(icat), Template{icat}, ploco(icat));
    end
end

% revisit the snippets and time align them w their median snippet
CatBorder = 1+cumsum([0 EventCount]);
NsamSpikeWidth = round(SpikeWidth/dt);
Wsn = exp(-0.5*(2*t0snip/SpikeWidth).^2);
icat = 0; ireject = [];
for ievent=1:Nevent,
    if ievent==CatBorder(icat+1),
        if icat==Ncat, break; end
        icat = icat+1
        while isequal(0, EventCount(icat)), icat=icat+1; end
        WmedSnip = Wsn.*Template{icat};  % center-windowed median snip of current category
        AllSnippets{icat} = [];
    end
    tsnip = tmax(ievent) + t0snip;
    snip = interp1(t,real(y),tsnip,'linear',0);
    Wysnip = Wsn.*snip;
    [XC, ilag] = xcorr(WmedSnip, Wysnip, NsamSpikeWidth, 'coeff');
    [MaxCor, imax] = max(XC);
    delta_t = dt*ilag(imax);
    %     if ~isequal(0,tshift),
    %         f6;
    %         plot(t0snip, WmedSnip,'k');
    %         xplot(t0snip, Wysnip, 'g');
    %         xplot(t0snip+delta_t, Wysnip, 'r'); pause
    %     end
    csnip(ievent,1) = MaxCor;
    if MaxCor>=CorCrit,
        tmax(ievent) = tmax(ievent) + delta_t;
        AllSnippets{icat} = [AllSnippets{icat}; snip];
    else,
        ireject = [ireject ievent]; % black list
    end
end
catIdx(ireject)=0;

% re-compute median snippet for each magnitude class
for icat=1:Ncat,
    if EventCount(icat)>0,
        Template{icat} = median(AllSnippets{icat},1);
    end
end

% sort events according to time order.
[tmax, ymax, catIdx, csnip] = sortAccord(tmax, ymax, catIdx, csnip,tmax);

% ID info; Use SGSR-style info if available
if isfield(D,'ExpID'),
    ExpID = D.ExpID; RecID = D.RecID; icond = D.icond; 
else,
    ExpID = ''; RecID = ''; icond = []; 
end
ABFname = D.Header.XX_FullFileName;
SweepDur = D.sweepDur_ms;
% collect all relevant (?) info in outpt arg struct S
S = CollectInStruct(ABFname, ExpID, RecID, icond, dt, SweepDur, '-', ...
    isweep, ichan, Vmax, '-', ...
    EventCount, Ncat, SpikeWidth, Ymean, '-', ...
    Template, t0snip, tmax, ymax, catIdx, csnip, ireject);

Thr = AlignThr;
if isempty(Thr),
    Thr = S.Ymean+min(S.Vmax)/2;
end
S.tshift = templot(S,Thr,0); % last zero arg: don't plot
for icat=1:Ncat,
    tsh = S.tshift(icat);
    if ~isnan(tsh),
        iev = find(S.catIdx==icat);
        S.tmax(iev) = S.tmax(iev)+tsh; 
    end
end
S.AlignThr = Thr;
S = xspikeAlign(S,-0.5, -0.2);

putcache(CFN, 1e4, CacheParam, S);
xspikeplot(S);














