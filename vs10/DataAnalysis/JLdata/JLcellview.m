function JLcellview(igerbil, icell);
% JLcellview - cell viewer for JL data
%   JLcellview(igerbil, icell)
%   or
%   JLcellview(icell_run)
%   or 
%   JLcellview(Uidx)

DB = JLdbase;
if nargin==1 && igerbil<100, %   JLcellview(icell_run)
    icell_run = igerbil;
elseif nargin==1, %   JLcellview(Uidx)
    icell_run = getfield(JLdbase(igerbil), 'icell_run');
else, % JLcellview(igerbil, icell)
    idb = find([DB.iexp]==igerbil & [DB.icell]==icell ,1);
    if isempty(idb),
        error(['Experiment ' num2str(igerbil) 'does not have data for Cell # ' num2str(icell) '.']);
    end
    icell_run = DB(idb).icell_run;
end


CanBeUsed = local_xlsfield(icell_run,'N');
ShouldBeUsed = local_xlsfield(icell_run,'P');

[Hd, Str, Uidx] = local_listbox(DB, icell_run, CanBeUsed, ShouldBeUsed);
figh = local_GUI(Hd, Str, Uidx);

%==============================================
function T = local_sep(S, isep, sep);
isep = find(isep);
Nsep = numel(isep);
isep = [0, isep(:).', size(S,1)];
SEP = repmat(sep, 1, size(S,2));
T = [];
for ii=1:Nsep+1,
    Block = S(isep(ii)+1:isep(ii+1),:);
    T = [T; Block; SEP];
end
T(end,:) = [];

function figh = local_GUI(Hd, Str, Uidx);
db = JLdbase(Uidx(1)); % example element of JLdbase
Xmax = 900; Ymax = 600;
PushButton = struct('style', 'pushbutton', 'fontsize', 12);
AnaButton = @(figh, y, Str, CB)uicontrol(figh, PushButton, 'position', [Xmax-100 y 90 30], 'String', Str, 'callback', CB);
figPos = [30 127 Xmax Ymax];
figh = newGUI('Cellview', sprintf('exp %d  cell %d    (icell_run = %d)', db.iexp, db.icell, db.icell_run), {@JLcellview db.iexp db.icell});
set(figh ,'units', 'pixels', 'position', figPos, 'visible', 'on');
set(figh ,'units', 'normalized');
Header = uicontrol(figh, 'style' ,'text', 'position', [20 Ymax-120 Xmax-140 15], ...
    'fontname', 'courier' ,'fontsize', 10, 'string', Hd, 'HorizontalAlignment', 'left');
RecList = uicontrol(figh ,'style', 'listbox', 'position', [20 20 Xmax-130 Ymax-160], ...
    'fontname', 'courier' ,'fontsize', 10, 'callback' ,@local_listboxmove,...
    'keypressfcn', @local_keypress, 'string', Str);
uicontrol(figh, 'position', [Xmax-350 Ymax-110 345 110], local_cellQuality(db.icell_run));
%================================buttons=================================
Yup = Ymax-50; % Y pos of upper row of buttons
spontRMSButton = uicontrol(figh, PushButton, 'position', [80 Yup 100 30], ...
    'Callback', @local_spontRMS, 'string', 'spontRMS', 'tooltipstring', 'Plot magnitude of spont inputs over time.');
RecIndexButton = uicontrol(figh, PushButton, 'position', [200 Yup 100 30], ...
    'Callback', @local_recindex, 'string', 'RecordingIndex', 'tooltipstring', 'Display unique recording index in command window.');
RecycleButton = uicontrol(figh, PushButton, 'position', [320 Yup 100 30], ...
    'Callback', @local_recfield, 'string', 'RecField', 'tooltipstring', 'Receptive-field view of periodograms.');
PowerButton = uicontrol(figh, PushButton, 'position', [440 Yup 100 30], ...
    'Callback', @local_powerfield, 'string', 'PowerField', 'tooltipstring', 'Receptive-field view of power-law plots.');
%=== right column of analysis buttons
ViewRecButton = AnaButton(figh, Ymax-140, 'View rec', @local_ViewRec);
TwoDimButton = AnaButton(figh, Ymax-180, '2D', @local_2D);
CycleButton = AnaButton(figh, Ymax-220, 'Cycles', @local_cycle);
SpikeButton = AnaButton(figh, Ymax-260, 'Spikes', @local_spikes);
TripletButton = AnaButton(figh, Ymax-300, 'Triplet', @local_triplet);
PhaseButton = AnaButton(figh, Ymax-340, 'Phase/freq', @local_phase_freq);
APthrButton = AnaButton(figh, Ymax-380, 'APthr', @local_APthr);
PeaksButton = AnaButton(figh, Ymax-420, 'peaKs', @local_peaks);
PeaksButton = AnaButton(figh, Ymax-460, 'AP Latency', @local_lat);
DYb = 55; Ygap = 260;
%==pulldown menus====
Data = uimenu(figh, 'label', '&Data');
UpdateItem = uimenu(Data, 'label', '&Update', 'Accelerator', 'U', 'callback', @local_update);
RecIndexItem = uimenu(Data, 'label', '&Paramview', 'Accelerator', 'P', 'callback', @local_paramview);
ListItem = uimenu(Data, 'label', '&Display List', 'Accelerator', 'D', 'callback', @local_displaylist);
Ana = uimenu(figh, 'label', '&Analysis');
ViewRecItem = uimenu(Ana, 'label', '&View rec', 'Accelerator', 'V', 'callback', @local_recindex);
TwoDimItem = uimenu(Ana, 'label', '2&D', 'Accelerator', 'D', 'callback', @local_2D);
CycleItem = uimenu(Ana, 'label', '&Cycles', 'Accelerator', 'C', 'callback', @local_cycle);
SpikesItem = uimenu(Ana, 'label', '&Spikes', 'Accelerator', 'S', 'callback', @local_spikes);
TripletItem = uimenu(Ana, 'label', '&Triplet', 'Accelerator', 'T', 'callback', @local_triplet);
PhaseItem = uimenu(Ana, 'label', '&Phase/freq', 'Accelerator', 'P', 'callback', @local_phase_freq);
APthrItem = uimenu(Ana, 'label', '&Apthr', 'Accelerator', 'A', 'callback', @local_APthr);
PeaksthrItem = uimenu(Ana, 'label', 'pea&Ks', 'Accelerator', 'K', 'callback', @local_peaks);
handles = CollectInStruct(figh, RecList, RecIndexButton, ...
    ViewRecButton, TwoDimButton, CycleButton, SpikeButton, '-', ...
    Data, UpdateItem, RecIndexItem, ListItem, '-',...
    Ana, ViewRecItem, TwoDimItem, CycleItem, TripletItem, PhaseItem, SpikesItem);
setGUIdata(figh, 'handles', handles); % remember the handles in figure, so callbacks may use them
setGUIdata(figh, 'Uidx', Uidx); % remember the unique rec indices
setGUIdata(figh, 'db', db); % remember the example database element

% =======helpers====================
function [figh, varargout] = local_handle(h, varargin);
% get named handle starting from any handle in GUI
figh = parentfigh(h); % GUI fig handle
H = getGUIdata(figh, 'handles');
for ii=1:numel(varargin),
    varargout{ii} = H.(varargin{ii});
end

function uidx = local_uidx(Src, Ev);
% get nique rec index of currently selected rec of listbox
uidx = getGUIdata(parentfigh(Src), 'Uidx'); % rec indices of whole list
[figh, hr] = local_handle(Src, 'RecList');
uidx = uidx(get(hr, 'value'));


%==============actual analysis fncs===========================
function local_spontRMS(Src, Ev);
figh = parentfigh(Src);
Uidx = denan(getGUIdata(figh, 'Uidx'));
D = JLdatastruct(Uidx);
Cst = JLcycleStats(Uidx);
Tim = [D.MinutesSinceRecStart];
DT = min(diff(Tim))/2;
[D, Cst] = sortAccord(D, Cst, Tim);
inewstim = find(diff([D.iseries 1e7])~=0);
set(figure,'units', 'normalized', 'position', [0.19 0.608 0.754 0.291]);
plot(Tim, sqrt([Cst.VarSpont1]), 'o-b');
xplot(Tim, sqrt([Cst.VarSpont2]), '*-b');
xplot(Tim, sqrt([Cst.VarDriftSpont1]), 'o-r');
xplot(Tim, sqrt([Cst.VarDriftSpont2]), '*-r');
ylim([0 1.2*max(ylim)]);
legend({['spont 1 (' num2str(round(Cst(1).Spont1Dur)) ' ms)'], ...
    ['spont 2 (' num2str(round(Cst(1).Spont2Dur)) ' ms)'] , ...
    'slow spont 1' , ...
    'slow spont 2' ...
    });
xlabel('Time since start (minutes)');
ylabel('RMS (mV)');
title(['exp ' num2str(D(1).iexp) ' cell ' num2str(D(1).icell)]);
% mark series chunks
fenceplot([Tim(inewstim)-DT], ylim, 'color', 'k');
for ii=1:numel(inewstim),
    id = inewstim(ii);
    iseries = D(id).iseries;
    tmean = mean([D([D.iseries]==iseries).MinutesSinceRecStart]);
    text(tmean, 0.8*max(ylim), sprintf('%s %d dB', D(id).chan, D(id).SPL), ...
        'Rotation', 90, 'color', [0 0.6 0]);
end
% mark currently selected rec
uidx = local_uidx(Src);
if ~isnan(uidx),
    D = JLdatastruct(uidx);
    fenceplot(D.MinutesSinceRecStart, ylim, 'color', [0 0.7 0], 'linewidth',3);
end

function local_ViewRec(Src, Ev);
uidx = local_uidx(Src);
figure;
JLviewrec(uidx);

function local_2D(Src, Ev);
uidx = local_uidx(Src);
JLanova2(uidx);

function local_cycle(Src, Ev);
uidx = local_uidx(Src);
JLbeatVar2(uidx);

function local_spikes(Src, Ev);
uidx = local_uidx(Src);
JLspikes(uidx);

function local_triplet(Src, Ev);
uidx = local_uidx(Src);
try, 
    JLbinint2(uidx);
catch,
    errordlg(lasterr, 'No triplets for you')
end

function local_phase_freq(Src, Ev);
uidx = local_uidx(Src);
JLpeakPhase(-uidx, 5.5);


function local_recfield(Src, Ev);
uidx = local_uidx(Src);
JLrecyclefield(uidx);

function local_powerfield(Src, Ev);
uidx = local_uidx(Src);
JLpowerfield(uidx);

function local_recindex(Src, Ev);
uidx = local_uidx(Src)

function local_APthr(Src, Ev);
uidx = local_uidx(Src);
JLspikeThr(uidx);

function local_peaks(Src, Ev);
uidx = local_uidx(Src);
JLinterpeak(uidx);

function local_lat(Src, Ev);
uidx = local_uidx(Src);
qq = JL_APana(uidx,1); % 1: do plot
%EPSPsize = qq.maxEslope; 

%=================================================================
function local_keypress(Src, Ev);
%Ev
figh = parentfigh(Src);
switch Ev.Key,
    case 'escape';
        commandwindow;
    case 'delete';
        close(figh);
    case '2'; % 
        local_2D(Src, Ev);
    case '3'; % 
        local_custom_3(figh, Ev);
    case '4'; % 
        local_custom_4(figh, Ev);
    case '5'; % 
        local_custom_5(figh, Ev);
    case '6'; % 
        local_custom_6(figh, Ev);
    case 'd'; % 
        local_2D(Src, Ev);
    case 'v'; % 
        local_ViewRec(figh, Ev);
    case 'c'; % 
        local_cycle(figh, Ev);
    case 'u'; % 
        local_update(figh, Ev);
    case 's'; % 
        local_typelog(figh, Ev);
    case 'p'; % 
        local_phase_freq(figh, Ev);
    case 'a'; % 
        local_APthr(figh, Ev);
end
drawnow;


function [CLR TC] = local_cell_color(icell_run);
TC = [0 0 0]; % black text
switch icell_run
    case {1 4 6 7 13 14 15 16 17 18 19 20 22 23 25 30 31 32 33 35 37 38 40 41 42},
        CLR = [0 1 0];
    case {2 3 5 9 10 26 29},
        CLR = [0 0 1];
        TC = [1 1 1];
    case {8 21 24 34 36 39 43 44 45},
        CLR = [1 0 0];
    case {11 12 28},
        CLR = [1 0.7 0.7];
end


function T = local_tooltip(ic);
[qq, qqt] = xlsread(FFN,1, Xr);
T = ['Usable: ' qqt{1}];

function S = local_cellQuality(ic);
UseTxt = local_xlsfield(ic, 'M', 'Useful');
StabTxt = local_xlsfield(ic, 'J', 'Stability');
BaseTxt = local_xlsfield(ic, 'K', 'Baseline var');
BinTxt = local_xlsfield(ic, 'R', 'Binaurality');
DrivTxt = local_xlsfield(ic, 'Z', 'Driven');
[BC, TC] = local_cell_color(ic);
S = struct('style', 'text', ...
    'fontsize', 10, ...
    'horizontalalign', 'left', ...
    'foregroundcolor', TC, ...
    'backgroundcolor', BC, ...
    'string', [UseTxt char(10) StabTxt char(10) BaseTxt char(10) BinTxt char(10) DrivTxt char(10)]);

function T = local_xlsfield(ic,C, Pref);
FFN = 'D:\USR\Gerard\MSO experiments list 20110923 GB.xls';
Xr = @(C)[C num2str(ic+1) ':' C num2str(ic+1)];
[dum, T] = xlsread(FFN,1, Xr(C));
if isempty(T), 
    T = '-';
else,
    T = T{1};
end
if nargin>2,
    T = [Pref ': ' T];
end

function [DB, cpref] = local_restrict(DB, F, Prf);
N = numel(DB);
if isequal('-', Prf),
    ipref = [];
else,
    Prf = strrep(Prf,' ','');
    Prf = strrep(Prf,'-',':');
    ipref = eval(Prf);
end
cpref = [blanks(N)' blanks(N)'];
cpref(ipref,1) = '+';
switch lower(F),
    case {'all' '-'}; % no restriction
    otherwise,
        F = strrep(F,' ','');
        F = Words2cell(F,',');
        isel = [];
        for ii=1:numel(F),
            f = F{ii};
            f = strrep(f,'-',':');
            f = strrep(f,'end',num2str(N));
            isel = [isel eval(f)];
        end
        DB = DB(isel);
        cpref = cpref(isel,:);
end

function [Hd, Str, Uidx] = local_listbox(DB, icell_run, CanBeUsed, ShouldBeUsed);
% restrict database to this cell and select relevant fields
DB = DB([DB.icell_run]==icell_run);
DB = sortAccord(DB, [DB.MinutesSinceRecStart]); % sort according to abs recording time
[DB.irec] = DealElements(1:numel(DB)); % add simple counting index
DB = DB([DB.freq] ~= 50);
[DB, cpref] = local_restrict(DB, CanBeUsed, ShouldBeUsed);
Sdisp = structpart(DB, {'irec' 'chan', 'freq', 'SPL'});
[Sdisp.Minutes] = DealElements(round(10*[DB.MinutesSinceRecStart])/10);
[Sdisp.series] = DealElements([DB.iseries]);
[Sdisp.series_run] = DealElements([DB.iseries_run]);
% store uniq rec indices
Uidx = [DB.UniqueRecordingIndex];
%structview(Sdisp);

% construct char matrix for listbox view
[Str, Hd] = struct2char(Sdisp, 8);
Str = [cpref Str];
%Str = [blanks(size(Str,1)).', Str]; % prepend spaces to disable navigation by key presses
Hd = [blanks(size(Hd,1)).', Hd];
% insert separators in display string & in rec indices, so they stay matched 
qsep = diff([Sdisp.series])~=0; % where to insert 
Str = local_sep(Str, qsep, '-');
Uidx = local_sep(Uidx(:), qsep, nan).';



function local_listboxmove(varargin);
% dummy