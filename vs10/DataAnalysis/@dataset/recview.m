function H = recview(DS, Chan, icond, irep);
% Dataset/recview - view single traces of analog data
%    recview(DS, Chan) plots single traces of analog data from channel Chan
%    of dataset DS. Chan may be a number or the full name (Dataset/anachan).
%
%    recview(DS, Chan, icond, irep) selects stimulus condition icond and
%    repetition irep. Defaults are Chan=1, icond=1, irep=1.
%
%    See also Dataset/anaview, Dataset/anadata, Dataset/anachan.
global dc_offset;
[Chan, icond, irep] = arginDefaults('Chan,icond,irep', 1,1,1);

if isa(DS, 'dataset'), % initialize plot GUI
    local_init(DS, Chan, icond, irep);
end





%=================================
function he = local_query(figh, Pos, DX, Prompt, Val, CB, TT);
CLR = get(figh, 'color');
uicontrol(figh, 'units', 'normalized', 'position', Pos, ...
    'backgroundcolor', CLR, 'style', 'text', 'callback', 'dir', 'fontsize', 10, ...
    'horizontalalign', 'right', 'string', Prompt, 'tooltipstring', TT);
he = uicontrol(figh, 'units', 'normalized', 'position', Pos+[DX 0 0 0], ...
    'style', 'edit', 'callback', 'dir', 'fontsize', 10, 'horizontalalign', 'left', ...
    'backgroundcolor', 'w', 'string', num2str(Val), 'callback', CB);

function local_toggle(Src, Ev);
Str = get(Src, 'userdata');
cstr = get(Src, 'string');
N = numel(Str);
ii = strmatch(cstr, Str);
ii = 1+rem(ii,N);
set(Src, 'string', Str{ii});
local_plot(Src, Ev);

function hb = local_togglebutton(figh, Pos, Str);
hb = uicontrol(figh, 'style', 'pushbutton', 'units', 'normalized', ...
    'position', Pos, 'callback', @local_toggle, 'userdata', Str, 'String', Str{1}, ...
    'fontsize', 10);

function hb = local_actionbutton(figh, Pos, Str, CB);
hb = uicontrol(figh, 'style', 'pushbutton', 'units', 'normalized', ...
    'position', Pos, 'callback', CB, 'String', Str, 'fontsize', 10);


function local_init(DS, Chan, icond, irep);
figh = figure;
figh = figh.Number;
set(figh, 'units', 'normalized', 'position', [0.0156 0.0459 0.8461 0.5]);
ha = axes('Position', [0.0156 0.5 0.8461 0.35], 'fontsize', 10);
ha2 = axes('Position', [0.0156 0.1 0.8461 0.35], 'fontsize', 10);
Ymax = 0.83; DY = 0.1;
H.paramview = local_actionbutton(figh, [0.6 0.93 0.08 0.07], 'stimulus', @(Src,Ev)paramview(DS));
H.Chan = local_query(figh, [0.86 Ymax-0*DY 0.05 0.07], 0.06, 'chan:', Chan, @local_plot, 'AD channel.');
H.icond = local_query(figh, [0.86 Ymax-1*DY 0.05 0.07], 0.06, 'icond:', icond, @local_plot, 'Stimulus condition index.');
H.irep = local_query(figh, [0.86 Ymax-2*DY 0.05 0.07], 0.06, 'irep:', irep, @local_plot, 'Repetition index.');
H.flt50Hz = local_togglebutton(figh, [0.88 Ymax-3*DY 0.1 0.07], {'no filter' 'reject 50 Hz'});
H.dt = local_query(figh, [0.86 Ymax-4*DY 0.05 0.07], 0.06, 'dt:', 0.15, @local_plot, 'Smoothing time (ms) for time derivative dV/dt.');
H.mode = local_togglebutton(figh, [0.92 Ymax-5*DY 0.05 0.07], {'V' 'dV/dt'});
H.events = local_togglebutton(figh, [0.89 Ymax-6*DY 0.09 0.07], {'no APs' 'stored APs' 'pick APs'});
H.APthr = local_query(figh, [0.86 Ymax-7*DY 0.05 0.07], 0.06, 'thr:', [], @local_plot, 'AP threshold (mV).');
H.paramview = local_actionbutton(figh, [0.86 Ymax-8*DY 0.08 0.07], 'peak stats', @local_peakstats);
H.ExportGraph = local_actionbutton(figh, [0.94 Ymax-8*DY 0.08 0.07], 'Export', @local_Export);

setGUIdata(figh, 'dataset', DS);
setGUIdata(figh, 'hPlotAxes', ha);
setGUIdata(figh, 'hStimAxes', ha2);
setGUIdata(figh, 'Viewparam', CollectInStruct(Chan, icond, irep));
setGUIdata(figh, 'EditHandles', H);
set(figh, 'toolbar', 'figure')
local_plot(H.Chan, nan, 1);


function X = local_read(H);
FNS = fieldnames(H);
for ii=1:numel(FNS),
    fn = FNS{ii};
    h = H.(fn);
    s = get(h,'string');
    hstyle = get(h, 'style');
    if isequal('edit', hstyle),
        X.(fn) = str2num(s);
    elseif isequal('pushbutton', hstyle),
        X.(fn) = s;
    end
end


function local_plot(Src, Ev, firsttime);
if nargin<3,firsttime=0; end
figh = parentfigh(Src);
figh = figh.Number;
H = getGUIdata(figh, 'EditHandles');
X = local_read(H);
ds = getGUIdata(figh, 'dataset');
ha = getGUIdata(figh, 'hPlotAxes');
axes(ha);
cla(ha);
ylim auto;
iplot = 0;
if X.irep == -1
    local_plot_average(ds,firsttime,Src);
    local_plot_stim(ds,firsttime,Src);
    return
end
for ic = X.icond(:).',
    for ir = X.irep(:).',
        iplot = iplot+1;
        [D dt t0] = anadata(ds, X.Chan, ic, ir);
        if isequal('reject 50 Hz', X.flt50Hz),
            D = reject50Hz(dt,D);
        end
        if isequal('pick APs', X.events),
            if ~isscalar(X.APthr)
                waitfor(warndlg('No AP threshold specified.','No input'));
                axes(ha);
                X.APthr = max(D);
            end
            SPT(ic,ir) = peakpicker([dt t0],D, X.APthr);
        end
        % raw trace or derivative?
        if isequal('dV/dt', X.mode),
            D = smoothen(D, X.dt, -dt);
            D = diff(D)/dt;
        end
        xdplot([dt t0], D, lico(iplot));
    end
end
% unlock & lock x and y limits
if firsttime, 
    xlim([0 dt*numel(D)]); 
    XL = xlim;
else,
    XL = xlim;
    xlim auto;
end
if isequal('pick APs', X.events) && isscalar(X.APthr),
    xplot(xlim, X.APthr*[1 1], 'k--');
end
ylim(ylim); % fix
xlim(XL);
xlabel('time (ms)', 'fontsize', 10);
ylabel('rec (V)', 'fontsize', 10);
CondStr = CondLabel(ds.Stim.Presentation, X.icond);
if iscell(CondStr), CondStr = cell2words(CondStr([1 end]), '...'); end;
title([IDstring(ds, 'full') ' ' CondStr], 'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');

TracePlotInterface(figh);
% add events?
if isequal('stored APs', X.events),
    SPT = spiketimes(ds,1, 'no-unwarp');
end
if ~isequal('no APs', X.events),
    iplot = 0;
    for ic = X.icond(:).',
        for ir = X.irep(:).',
            iplot = iplot+1;
            spt = SPT{ic,ir};
            LC = lico(iplot);
            LC.color = 0.5+LC.color*0.5;
            LC.linewidth = 2;
            LC.linestyle = '--';
            fenceplot(spt, ylim, LC);
        end
    end
end


function local_peakstats(Src, Ev);
figh = parentfigh(Src);
figh = figh.Number;
H = getGUIdata(figh, 'EditHandles');
X = local_read(H);
ds = getGUIdata(figh, 'dataset');
[D dt] = anadata(ds, X.Chan, X.icond, X.irep);
if isequal('reject 50 Hz', X.flt50Hz),
    D = reject50Hz(dt,D);
end
set(figure,'units', 'normalized', 'position', [0.688 0.171 0.264 0.264])
peakstats(dt, D, 0.5, 100);

function local_plot_stim(ds,firsttime,Src)
figh = parentfigh(Src);
figh = figh.Number;
ha = getGUIdata(figh, 'hStimAxes');
axes(ha);
cla;


H = getGUIdata(figh, 'EditHandles');
X = local_read(H);
icond = X.icond(1);
chan = X.Chan;

W = ds.Stim.Waveform(icond,chan);
dt = 1e3/W(1).Fsam; % sample period in ms
W = samples(W);
set(gca,'NextPlot','add')

% Plotting

x = W;
clr = get(0,'defaultAxesColorOrder');
xdplot([dt 0],real(x), 'color', clr(end,:));


xlim('auto');
ylim('auto');



function local_plot_average(ds,firsttime,Src)
figh = parentfigh(Src);
figh = figh.Number;
H = getGUIdata(figh, 'EditHandles');
X = local_read(H);
ds = getGUIdata(figh, 'dataset');
ha = getGUIdata(figh, 'hPlotAxes');
axes(ha);
cla(ha);
ylim auto;
iplot = 0;
for ic = X.icond(:).',
    D_avg = [];
    for ir = 1:ds.Stim.Presentation.Nrep,
        iplot = iplot+1;
        [D dt t0] = anadata(ds, X.Chan, ic, ir);
        D_avg = [D_avg D]; 
        if isequal('reject 50 Hz', X.flt50Hz),
            D = reject50Hz(dt,D);
        end
        if isequal('pick APs', X.events),
            if ~isscalar(X.APthr)
                waitfor(warndlg('No AP threshold specified.','No input'));
                axes(ha);
                X.APthr = max(D);
            end
            SPT(ic,ir) = peakpicker([dt t0],D, X.APthr);
        end
        % raw trace or derivative?
        if isequal('dV/dt', X.mode),
            D = smoothen(D, X.dt, -dt);
            D = diff(D)/dt;
        end
        xdplot([dt t0], D, lico(iplot));
    end
    D_avg = mean(D_avg,2);
    setGUIdata(figh, 'dc_offset',mean(D_avg));
    setGUIdata(figh, 'range',range(D_avg));
    xdplot([dt t0], D_avg, lico(iplot+1),'LineWidth',5);
end
% unlock & lock x and y limits
if firsttime, 
    xlim([0 dt*numel(D)]); 
    XL = xlim;
else,
    XL = xlim;
    xlim auto;
end
if isequal('pick APs', X.events) && isscalar(X.APthr),
    xplot(xlim, X.APthr*[1 1], 'k--');
end
ylim(ylim); % fix
xlim(XL);
xlabel('time (ms)', 'fontsize', 10);
ylabel('rec (V)', 'fontsize', 10);
CondStr = CondLabel(ds.Stim.Presentation, X.icond);
if iscell(CondStr), CondStr = cell2words(CondStr([1 end]), '...'); end;
title([IDstring(ds, 'full') ' ' CondStr], 'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');

% TracePlotInterface(figh);
% add events?
if isequal('stored APs', X.events),
    SPT = spiketimes(ds,1, 'no-unwarp');
end
if ~isequal('no APs', X.events),
    iplot = 0;
    for ic = X.icond(:).',
        for ir = X.irep(:).',
            iplot = iplot+1;
            spt = SPT{ic,ir};
            LC = lico(iplot);
            LC.color = 0.5+LC.color*0.5;
            LC.linewidth = 2;
            LC.linestyle = '--';
            fenceplot(spt, ylim, LC);
        end
    end
end

function local_Export(Src, Ev, firsttime);
if nargin<3,firsttime=0; end
figh = parentfigh(Src);
figh = figh.Number;
H = getGUIdata(figh, 'EditHandles');
X = local_read(H);
ds = getGUIdata(figh, 'dataset');
ha2 = getGUIdata(figh, 'hStimAxes'); % plot stim on the other axis (Jan 2018)
ha = getGUIdata(figh, 'hPlotAxes');
ylim auto;
iplot = 0;
if X.irep == -1
    ExportMean = 1;
else
    ExportMean = 0;
end
D_tot = [];
for ic = X.icond(:).',
    D_exp = [];
    if X.irep == -1
        X.irep = 1:ds.Stim.Presentation.Nrep;
    end
    for ir = X.irep(:).',
        iplot = iplot+1;
        [D dt t0] = anadata(ds, X.Chan, ic, ir);
        if isequal('reject 50 Hz', X.flt50Hz),
            D = reject50Hz(dt,D);
        end
        if isequal('pick APs', X.events),
            if ~isscalar(X.APthr)
                waitfor(warndlg('No AP threshold specified.','No input'));
                axes(ha);
                X.APthr = max(D);
            end
            SPT(ic,ir) = peakpicker([dt t0],D, X.APthr);
        end
        % raw trace or derivative?
        if isequal('dV/dt', X.mode),
            D = smoothen(D, X.dt, -dt);
            D = diff(D)/dt;
        end
        D_exp = [D_exp D];
        
    end
    D_avg = mean(D_exp,2);
    D_exp = [D_exp D_avg];
    D_tot = [D_tot D_exp];
end 

icond = X.icond(1);
chan = X.Chan;

W = ds.Stim.Waveform(icond,chan);
dt = 1e3/W(1).Fsam; % sample period in ms

% Added to fix struct error (Jan 2018)
if iscell(W.Samples),
    stim_construct = [];
    for ii = 1:length(W.Samples)
    stim_part = W.Samples{ii};
    stim_part_rep = W.Nrep(ii);
    stim_construct = [stim_construct;repmat(stim_part,stim_part_rep,1)];
    end
    W = stim_construct;
else
    W = W.Samples;
end
% W = samples(W);
set(gca,'NextPlot','add')

% Plotting

x = W;
clr = get(0,'defaultAxesColorOrder');
axes(ha2); % plot the stim on the second gui axis (Jan 2018)
xdplot([dt 0],real(x), 'color', clr(end,:));

Export.Time = (dt*(0:N-1));
Export.Recordings = D_tot;
Export.Stim = W;
Export.Fsam_rec = ds.Rec.RecordInstr(1).Fsam;
Export.Fsam_stim = ds.Stim.Fsam;
save(IDstring(ds),'Export');