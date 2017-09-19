function PlotSpikeMetrics(M, ievent);
% PlotSpikeMetrics - plot trace annotated with spike metrics 
%    PlotSpikeMetrics(X,M,ievent);

persistent hfig

if nargin<2, ievent=[]; end; % default: don't jump to an event

% does a spikemetric plot exist having the same M?
shhv=get(0,'ShowHiddenHandles'); set(0,'ShowHiddenHandles','on');
hfig = findobj('type', 'figure', 'tag', mfilename, 'userdata', M);
set(0,'ShowHiddenHandles',shhv);
if isempty(hfig), % no - create new figure
    hfig = figure('integerhandle', 'off', ...
        'numbertitle', 'off');
    set(hfig,'units', 'normalized', 'position', [0.0445 0.373 0.952 0.499]);
    local_initplot(hfig, M);
else, % yes - bring to the fore
    figure(hfig);
    set(hfig, 'handlevi', 'on'); % temporary access 
end

if ~isempty(ievent),
    local_jump2event(ievent);
end
set(hfig, 'handlevi', 'callback'); % end tem access

%==============================================================
function local_CallBack(Src,Ev,varargin);
Shifted = ismember('shift', Ev.Modifier);
switch Ev.Key,
    case 'leftarrow', % move to left
        if Shifted, Factor = 1; else, Factor = 0.1; end
        xlim(xlim-0.5*Factor*diff(xlim)); pause(0.025);
        xlim(xlim-0.5*Factor*diff(xlim)); 
    case 'rightarrow', % move to right
        if Shifted, Factor = 1; else, Factor = 0.1; end
        xlim(xlim+0.5*Factor*diff(xlim)); pause(0.025);
        xlim(xlim+0.5*Factor*diff(xlim)); 
    case 'uparrow',  % zoom in
        if Shifted, Factor = 0.5; else, Factor = 0.9; end
        xlim(mean(xlim)+0.5*sqrt(Factor)*[-1 1]*diff(xlim)); pause(0.025);
        xlim(mean(xlim)+0.5*sqrt(Factor)*[-1 1]*diff(xlim));
    case 'downarrow', % zoom out
        if Shifted, Factor = 0.5; else, Factor = 0.9; end
        xlim(mean(xlim)+0.5/sqrt(Factor)*[-1 1]*diff(xlim)); pause(0.025);
        xlim(mean(xlim)+0.5/sqrt(Factor)*[-1 1]*diff(xlim));
    case 'home', % to start of trace
        xlim(xlim-min(xlim));
    case 'end', % to end of trace
        DX = diff(xlim);
        xlim('auto');
        xlim(max(xlim)+[-DX 0]);
    case 'a', % default view
        xlim('auto');
    case '0',
        xlim([0 30]);
    case 'g',
        ievent = str2double(inputdlg('Specify event #'));
        local_jump2event(ievent);
    case 's',
        M=MyFlag(mfilename); % retrieve M input arg
        XL = xlim; DX=diff(XL); X0=XL(1);
        if ~Shifted, % jump to next onset
            X0 = min(M.onset(M.onset>=X0+1.1));
        else, % previous
            X0 = max(M.onset(M.onset<=X0-1.1));
        end
        if ~isempty(X0),
            xlim((X0-1)+DX*[0 1]);
        end
    otherwise, 
        %disp(Ev)
end

function  local_initplot(hfig, M);
% get complete trace
D = getABFfromSpikeStruct(M.X); % retrieve ABF recording
dt = M.X.dt;
Trace = cat(1,D.AD(M.X.isweep, M.X.ichan).samples); % concatenate the sweeps
T = timeaxis(Trace,dt);
% create figure
set(hfig, 'tag', mfilename, 'userdata', M);
plot(T, Trace, 'k', 'linewidth', 2); YL = ylim';
xlim([0 50]); drawnow;
fenceplot(M.tMaxEPSPrate, YL, 'b');
xplot(M.tMaxAP, M.MaxAP, 'r*');
xplot(M.tMinAP, M.MinAP, 'ro');
xplot([1;1]*[M.tMaxPrespike], [M.MaxPrespike; M.MaxPrespike-M.SizePrespike], 'm.-', 'markersize', 12);
ylim(min(ylim)+1.1*[0 1]*diff(ylim)); % fix Ylimits to prevent rescaling when browsing using keys
YL = ylim; DY = diff(YL);
dashplot(M.tMaxEPSP, M.MaxEPSP, 0.3, 'color', [0.7 0 0]);
dashplot(M.onset+0.5*M.burstDur, 0*M.onset+YL(1)+0.02*DY, M.burstDur, ...
    'color', [0.9 0.8 0.8], 'linewidth', 8);
for ie=1:M.Nev,
    yoffs = 0.04*rem(ie,2)*diff(ylim);
    text(M.tMaxEPSPrate(ie), min(ylim)+yoffs+0.9*diff(ylim), num2str(ie), ...
        'color', 'b', 'backgroundcolor', 'w', 'clipping', 'on', ...
        'horizontalalign', 'center');
end
xlabel('Time (ms)', 'fontsize', 12);
ylabel('V_{ex} (mV)', 'fontsize', 12);
set(hfig,'keypressfcn', {@local_CallBack});
IDstr = [M.ExpID '/' M.RecID '/' num2str(M.icond) ' ---- ' num2str(M.Xval) ' ' M.Xunit];
MyFlag(mfilename, M); % store M
title(IDstr);
set(hfig, 'name', ['Trace ' IDstr]);


function local_jump2event(ievent);
M = MyFlag(mfilename); % retrieve M input arg
if ~betwixt(ievent,0,M.Nev+1), return; end
DX = diff(xlim);
if DX>200, DX = 30; end
xlim(M.tMaxEPSPrate(ievent)+[-1 1]*DX/2);
yoffs = 0.04*rem(ievent,2)*diff(ylim);
text(M.tMaxEPSPrate(ievent), min(ylim)+yoffs+0.9*diff(ylim), num2str(ievent), ...
    'color', 'k', 'fontweight', 'bold', 'backgroundcolor', 'r', 'clipping', 'on', ...
    'horizontalalign', 'center');



