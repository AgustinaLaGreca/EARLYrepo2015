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
if ~isa(figh,'double')
    figh = figh.Number;
end;
set(figh, 'units', 'normalized', 'position', [0.0156 0.0459 0.8461 0.5]);
ha = axes('Position', [0.0156 0.5 0.8461 0.35], 'fontsize', 10);
ha2 = axes('Position', [0.0156 0.1 0.8461 0.35], 'fontsize', 10);
@ -92,6 +95,9 @@ end
function local_plot(Src, Ev, firsttime);
if nargin<3,firsttime=0; end
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
H = getGUIdata(figh, 'EditHandles');
X = local_read(H);
ds = getGUIdata(figh, 'dataset');
@ -170,6 +176,9 @@ end

function local_peakstats(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
H = getGUIdata(figh, 'EditHandles');
X = local_read(H);
ds = getGUIdata(figh, 'dataset');
@ -182,6 +191,9 @@ peakstats(dt, D, 0.5, 100);

function local_plot_stim(ds,firsttime,Src)
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
ha = getGUIdata(figh, 'hStimAxes');
axes(ha);
cla;
@ -211,6 +223,9 @@ ylim('auto');

function local_plot_average(ds,firsttime,Src)
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
H = getGUIdata(figh, 'EditHandles');
X = local_read(H);
ds = getGUIdata(figh, 'dataset');
@ -290,6 +305,9 @@ end
function local_Export(Src, Ev, firsttime);
if nargin<3,firsttime=0; end
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
H = getGUIdata(figh, 'EditHandles');
X = local_read(H);
ds = getGUIdata(figh, 'dataset');
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
W = samples(W);
set(gca,'NextPlot','add')

% Plotting

x = W;
clr = get(0,'defaultAxesColorOrder');
xdplot([dt 0],real(x), 'color', clr(end,:));

Export.Time = (dt*(0:N-1)).'
Export.Recordings = D_tot;
Export.Stim = W;
Export.Fsam_rec = ds.Rec.RecordInstr(1).Fsam;
Export.Fsam_stim = ds.Stim.Fsam;
save(IDstring(ds),'Export');
