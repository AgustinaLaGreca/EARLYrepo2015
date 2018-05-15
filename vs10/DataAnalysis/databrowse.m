function databrowse(Exp, varargin);
% databrowse - listing and standard analysis of recordings
%  databrowse('RG11308') opens a GUI which lists the stimulus conditions of
%  the recordings of experiment RG11308 and offers standard analyses.
%  varargin{1} is a cell string with the function and arguments which should
%  be called if databrowse is closed


if ischar(Exp),
    Exp = find(experiment(),Exp);
end
figh = local_GUI(Exp,varargin);
local_update(figh); % update listing of recordings


%==========================================
function figh = local_GUI(Exp, varargin);

%====DEFINITION=====
Xmax = 1200; Ymax = 600;
PushButton = struct('style', 'pushbutton', 'fontsize', 10);
AnaButton = @(figh, y, Str, CB)uicontrol(figh, PushButton, 'position', [Xmax-110 y 100 30], 'String', Str, 'callback', CB);
figPos = [30 127 Xmax Ymax];
figh = newGUI('Databrowse', name(Exp), {@databrowse name(Exp)});
if ~isempty(varargin);
    if ~ isempty(varargin{1})
        fcn = varargin{1};
        fcn = fcn{1};
        set(figh,'CloseRequestFcn',fcn);
    end
end
set(figh ,'units', 'pixels', 'position', figPos, 'visible', 'on');
set(figh ,'units', 'normalized');
RecList = uicontrol(figh ,'style', 'listbox', 'position', [10 20 Xmax-130 Ymax-100], ...
    'fontname', 'courier' ,'fontsize', 10, 'callback' ,@local_listboxmove,...
    'keypressfcn', @local_keypress);

%====BUTTONS=====
Yup = Ymax-50; % Y pos of upper row of buttons
ChanButton = uicontrol(figh, PushButton, 'position', [Xmax-100 Yup 60 30], ...
    'tooltipstring', 'AD channel specification of analog data. Click to toggle.');
ChanButton = double(ChanButton);
betoggle(ChanButton, {'AD-1', 'AD-2'  'AD-3'});
TypeLogButton = uicontrol(figh, PushButton, 'position', [80 Yup 100 30], ...
    'Callback', @local_typelog, 'string', 'typelog', 'tooltipstring', 'Display experiment log.');
TypeLogButton = double(TypeLogButton);
ParamViewButton = uicontrol(figh, PushButton, 'position', [200 Yup 100 30], ...
    'Callback', @local_paramview, 'string', 'paramview', 'tooltipstring', 'View stimulus parameters of recording (Ctrl-P).');
ParamViewButton = double(ParamViewButton);
UpdateButton = uicontrol(figh, PushButton, 'position', [320 Yup 100 30], ...
    'Callback', @local_update, 'string', 'update', 'tooltipstring', 'Update listing of recordings (Ctrl-U).');
UpdateButton = double(UpdateButton);
DYb = 40; % = right column of analysis buttons
DotRasterButton = AnaButton(figh, Ymax-3*DYb, 'dotraster', @local_dotraster);
DotRasterButton = double(DotRasterButton);
PowerSpecButton = AnaButton(figh, Ymax-4*DYb, 'power spectrum', @local_powerspec);
PowerSpecButton = double(PowerSpecButton);
ReflectanceButton = AnaButton(figh, Ymax-5*DYb, 'reflectance', @local_reflectance);
ReflectanceButton = double(ReflectanceButton);
ViewRecButton = AnaButton(figh, Ymax-6*DYb, 'view recording', @local_viewrec);
ViewRecButton = double(ViewRecButton);
DYb = 55; Ygap = 310;

%====CUSTOM BUTTONS===== % 6th button added by Jan 2018
for ibut=1:6
    eval(['CustomAna_' num2str(ibut) '_Button = AnaButton(figh, Ymax-Ygap-' num2str(ibut-1) '*DYb, ''...'', @local_dotraster);']);
end
% CustomAna_1_Button = AnaButton(figh, Ymax-Ygap, '...', @local_dotraster);
% CustomAna_2_Button = AnaButton(figh, Ymax-Ygap-DYb, '...', @local_dotraster);
% CustomAna_3_Button = AnaButton(figh, Ymax-Ygap-2*DYb, '...', @local_dotraster);
% CustomAna_4_Button = AnaButton(figh, Ymax-Ygap-3*DYb, '...', @local_dotraster);
% CustomAna_5_Button = AnaButton(figh, Ymax-Ygap-4*DYb, '...', @local_dotraster);

%==PULLDOWN MENUS====
Data = uimenu(figh, 'label', '&Data');
UpdateItem = uimenu(Data, 'label', '&update', 'Accelerator', 'U', 'callback', @local_update);
ParamViewItem = uimenu(Data, 'label', '&paramview', 'Accelerator', 'P', 'callback', @local_paramview);
ListItem = uimenu(Data, 'label', '&display list', 'Accelerator', 'D', 'callback', @local_displaylist);
Ana = uimenu(figh, 'label', '&Analysis');
RasterItem = uimenu(Ana, 'label', '&dotraster', 'Accelerator', 'R', 'callback', @local_dotraster);
PowerSpecItem = uimenu(Ana, 'label', 'power spectrum', 'Accelerator', 'S', 'callback', @local_powerspec);
ReflexItem = uimenu(Ana, 'label', 'reflectance', 'Accelerator', 'F', 'callback', @local_reflectance);
Custom_1_Item = uimenu(Ana, 'label', '&custom_1', 'Accelerator', '1', 'callback', @local_custom_1);
Custom_2_Item = uimenu(Ana, 'label', '&custom_2', 'Accelerator', '2', 'callback', @local_custom_2);
Custom_3_Item = uimenu(Ana, 'label', '&custom_3', 'Accelerator', '3', 'callback', @local_custom_3);
Custom_4_Item = uimenu(Ana, 'label', '&custom_4', 'Accelerator', '4', 'callback', @local_custom_4);
Custom_5_Item = uimenu(Ana, 'label', '&custom_5', 'Accelerator', '5', 'callback', @local_custom_5);
Export = uimenu(figh, 'label', '&Export');
ExportTable = uimenu(Export, 'label', '&DataBrowse Table', 'callback', @local_ExportTable);
ExportDotraster = uimenu(Export, 'label', '&DotRaster', 'callback', @(Src,Evt)local_ExportDotraster(Src,Evt,figh));
ExportSpectrum = uimenu(Export, 'label', '&Power Spectrum', 'callback', @(Src,Evt)local_ExportSpectrum(Src,Evt,figh));
ExportThresholdCurve = uimenu(Export, 'label', '&Threshold Curve', 'callback', @(Src,Evt)local_ExportThresholdCurve(Src,Evt,figh));
ExportPSTH = uimenu(Export, 'label', '&PSTH', 'callback', @(Src,Evt)local_ExportPSTH(Src,Evt,figh));
ExportCycleHisto = uimenu(Export, 'label', '&CycleHistogram', 'callback', @(Src,Ev)local_ExportCycleHisto(Src,Ev,figh));
ExportRatePlot = uimenu(Export, 'label', '&RatePlot', 'callback', @(Src,Ev)local_ExportRatePlot(Src,Ev,figh));
ExportResponsArea = uimenu(Export, 'label', '&ResponseArea', 'callback', @(Src,Ev)local_ExportResponsArea(Src,Ev,figh));
ExportFOISI = uimenu(Export, 'label', '&FOISI', 'callback', @(Src,Ev)local_ExportFOISI(Src,Ev,figh));
ExportCoeffvar = uimenu(Export, 'label', '&Coeffvar', 'callback', @(Src,Ev)local_ExportCoeffvar(Src,Ev,figh));

handles = CollectInStruct(figh, RecList, ChanButton, ParamViewButton, ...
    DotRasterButton, CustomAna_1_Button, CustomAna_2_Button, ...
    CustomAna_3_Button, CustomAna_4_Button, ...
    CustomAna_5_Button,CustomAna_6_Button, ... ,   '-', ...
    Data, UpdateItem, ParamViewItem, ListItem, '-',...
    Ana, RasterItem, Custom_1_Item, Custom_2_Item, '-',...
    Export, ExportTable, ExportDotraster, ExportSpectrum, ExportThresholdCurve, ...
    ExportPSTH, ExportCycleHisto, ExportRatePlot, ExportResponsArea, ...
    ExportFOISI, ExportCoeffvar); % 6th button added by Jan 2018
setGUIdata(figh, 'handles', handles); % remember the handles in figure, so callbacks may use them
setGUIdata(figh, 'Experiment', Exp); % remember the experiment in figure, so callbacks may use it
setGUIdata(figh, 'datagetter', getds(name(Exp))); % quick dataset retriever

function [figh, varargout] = local_handle(h, varargin);
% get named handle starting from any handle in GUI
figh = parentfigh(h); % GUI fig handle
if ~isa(figh,'double')
    figh = figh.Number;
end;
H = getGUIdata(figh, 'handles');
for ii=1:numel(varargin),
    varargout{ii} = H.(varargin{ii});
end

function local_ExportTable(Src,Ev)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for DataBrowse Table');
if isnumeric(FolderName)
    if FolderName == 0
        return;
    end
end
Experiment = current(experiment);
Table = ExportDataBrowseTable(Experiment);
save([FolderName filesep 'DataBrowseTable'],'Table');

function local_ExportDotraster(Src,Ev,figh)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for dotrqster data');
if isnumeric(FolderName)
    if FolderName == 0
        return;
    end
end
[ds, ichan] = local_curds(figh);
data_Struct = dotraster(ds);
save([FolderName filesep 'Dotraster_' IDstring(ds)],'data_Struct');

function local_ExportSpectrum(Src,Ev,figh)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for Spectrum data');
if isnumeric(FolderName)
    if FolderName == 0
        return;
    end
end
[ds, ichan] = local_curds(figh);
icond = local_icond(ds.Ncond);
if isempty(icond), return; end
data_struct = powerspec(ds, ichan, icond);
save([FolderName filesep 'PowerSpectrum_Cond' num2str(icond) '_' IDstring(ds)],'data_struct');

function local_ExportThresholdCurve(Src,Ev,figh)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for THR data');
[ds, ichan] = local_curds(figh);
data_struct = threshold_curve(ds);
save([FolderName filesep 'THRCurve_' IDstring(ds)],'data_struct');

function local_ExportPSTH(Src,Ev,figh)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for PSTH data');
if isnumeric(FolderName)
    if FolderName == 0
        return;
    end
end
[ds, ichan] = local_curds(figh);
data_struct = PSTH(ds);
save([FolderName filesep 'PSTH_' IDstring(ds)],'data_struct');

function local_ExportCycleHisto(Src,Ev,figh)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for CycleHistogram data');
if isnumeric(FolderName)
    if FolderName == 0
        return;
    end
end
[ds, ichan] = local_curds(figh);
data_struct = cyclehisto(ds);
save([FolderName filesep 'CycleHistogram_' IDstring(ds)],'data_struct');

function local_ExportRatePlot(Src,Ev,figh)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for RatePlot data');
if isnumeric(FolderName)
    if FolderName == 0
        return;
    end
end
[ds, ichan] = local_curds(figh);
data_struct = rateplot(ds);
save([FolderName filesep 'RatePlot_' IDstring(ds)],'data_struct');

function local_ExportResponsArea(Src,Ev,figh)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for ResponsArea data');
if isnumeric(FolderName)
    if FolderName == 0
        return;
    end
end
[ds, ichan] = local_curds(figh);
data_struct = resparea(ds);
save([FolderName filesep 'ResponsArea_' IDstring(ds)],'data_struct');

function local_ExportFOISI(Src,Evt,figh)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for FOISI data');
if isnumeric(FolderName)
    if FolderName == 0
        return;
    end
end
[ds, ichan] = local_curds(figh);
data_struct = FOISI(ds);
save([FolderName filesep 'FOISI_' IDstring(ds)],'data_struct');

function local_ExportCoeffvar(Src,Evt,figh)
FolderName = uigetdir('C:\ExpData\Exp','Select Destination folder for coeffvar data');
if isnumeric(FolderName)
    if FolderName == 0
        return;
    end
end
[ds, ichan] = local_curds(figh);
data_struct = coeffvar(ds);
save([FolderName filesep 'coeffvar_' IDstring(ds)],'data_struct');

function local_update(Src,Ev)
% update stimulus list and display it in the listbox
[figh, hreclist] = local_handle(Src, 'RecList');
oldNitem = size(get(hreclist ,'String'),1);
Exp = getGUIdata(figh, 'Experiment');
[LL RR TT] = stimlist(Exp);
Nrec = numel(TT);
setGUIdata(figh,'StimTypes', TT); % store stimulus types in GUI
RecList = [repmat(' ',[Nrec 1]) struct2char(LL) struct2char(RR)]; % heading space prevent side effects ...
Nitem = size(RecList,1);  %                                      ...  of pressing '1' etc. See local_keypress
if Nitem>oldNitem,
    lbt = get(hreclist, 'ListboxTop'); 
    lbv = get(hreclist, 'Value');
    lbv = max(lbv,1);
    newtop = max(Nitem-32, lbt);
    set(hreclist, 'string', RecList, 'ListboxTop', newtop);
    if newtop>lbv, set(hreclist, 'value', newtop); end
    local_listboxmove(figh, []);
end

function [ds, ichan] = local_curds(figh);
% Return the selected dataset and corresponding analog channel.
Exp = getGUIdata(figh, 'Experiment');
[figh, hreclist] = local_handle(figh, 'RecList');
irec = get(hreclist, 'Value');
DG = getGUIdata(figh, 'datagetter');
ds = DG(irec);
ichan = local_ADchan(figh, ds);
curds(ds); % "publish" ds

function ichan = local_ADchan(figh, ds);
% Return AD channel of the given dataset ds. 
% If only ADC chan one is present, do that one. Else use pushbutton.
dd = ds.Data;
has1 = isfield(dd, 'RX6_analog_1');
has2 = isfield(dd, 'RX6_analog_2');
if has1 && has2, % consult button
    [figh, hchanbut] = local_handle(figh, 'ChanButton');
    Str = get(hchanbut, 'String');
    ichan = str2num(Str(4:end));
elseif has1, ichan = 1;
elseif has2, ichan = 2;
else
%      error('No analog data in this dataset.'); % commented by Jan 2018
    ichan = nan;
end

function local_typelog(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
Exp = getGUIdata(figh, 'Experiment');
typelog(Exp);

function local_paramview(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
ds = local_curds(figh);
paramview(ds);

function local_displaylist(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[figh, hreclist] = local_handle(figh, 'RecList');
Str = get(hreclist, 'String');
more off;
disp(Str);

%======ANALYSIS============

% EVENTS
function local_dotraster(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
dotraster(ds);

function local_cyclehisto(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
cyclehisto(ds,figure);

function local_rateplot(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
rateplot(ds);

function local_resparea(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
resparea(ds);

function local_maskcorr(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
maskcorr(ds);

function local_PSTH(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
PSTH(ds);

function local_FOISI(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
FOISI(ds);

function local_AOISI(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
AOISI(ds);

function local_coeffvar(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
coeffvar(ds);

function local_revcor(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
revcor(ds);

% ANALOG DATA
function icond = local_icond(Ncond)
if Ncond==1, icond=1;
else, icond = inputdlg('iCond: '); 
    if ~isempty(icond),
        icond = str2num(icond{1});
    end
end

function local_anamean(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
IDstr = IDstring(ds, 'full');
icond = local_icond(ds.Ncond);
if ~icond  
    for icond=1:ds.Ncond
        [D dt t0] = anadata(ds,ichan,icond);
        D = anamean(ds, ichan, icond);
        [ScalFac, Unit] = conversionfactor(anachan(ds,ichan));
        D = D*ScalFac; % rec in true units (not volts of the ADC)
        xdplot([dt t0], D)
    end
else
    [D dt t0] = anadata(ds,ichan,icond);
    D = anamean(ds, ichan, icond);
    [ScalFac, Unit] = conversionfactor(anachan(ds,ichan));
    D = D*ScalFac; % rec in true units (not volts of the ADC)
    dplot([dt t0], D)
end
xlabel('time (ms)','fontsize',10);
ylabel(['rec (' Unit ')'],'fontsize',10);
title([IDstr '  AD-'  num2str(ichan) ' (' dataType(anachan(ds,ichan)) ')'],...
    'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');

function local_rmsplot(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
figure;
RMSplot(ds, ichan);

function local_supspec(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
icond = local_icond(ds.Ncond);
if isempty(icond), return; end
supspec(ds, ichan, icond);

function local_magn_phi(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
figure; Magn_Phase_Plot(ds, ichan);

function local_apple(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
set(figure,'units', 'normalized', 'position', [0.333 0.433 0.479 0.469])
apple(ds, ichan);

function local_powerspec(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
icond = local_icond(ds.Ncond);
if isempty(icond), return; end
figure; powerspec(ds, ichan, icond);

function local_reflectance(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
if ~isfield('PC_COM1',ds.Data), 
    errordlg('No PC_COM1 connected for this Experiment.');
else
    figure; plot(ds.Data.PC_COM1); title(IDstring(D, 'full'), 'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');
end

function local_viewrec(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
icond = local_icond(ds.Ncond);
if isempty(icond), return; end
recview(ds, ichan, icond);

function local_IOfun(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
ifreq = inputdlg('iFreq: '); 
if isempty(ifreq), return; end
ifreq = str2num(ifreq{1});
figure; IOfun(ds, ichan, ifreq);

function local_threshold_curve(Src, Ev);
figh = parentfigh(Src);
if ~isa(figh,'double')
    figh = figh.Number;
end;
[ds, ichan] = local_curds(figh);
threshold_curve(ds);

%LLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
function local_listboxmove(Src, Ev)
[figh, hlb, h1, h2, h3, h4, h5, h6] = local_handle(Src, 'RecList', ...
    'CustomAna_1_Button', 'CustomAna_2_Button', ...
    'CustomAna_3_Button', 'CustomAna_4_Button', ...
    'CustomAna_5_Button','CustomAna_6_Button');

irec = get(hlb, 'value');
StimTypes = getGUIdata(figh, 'StimTypes');
set(h1, 'String', '  ... ', 'callback', @nope);
set(h2, 'String', '  ... ', 'callback', @nope);
set(h3, 'String', '  ... ', 'callback', @nope);
set(h4, 'String', '  ... ', 'callback', @nope);
set(h5, 'String', '  ... ', 'callback', @nope);
set(h6, 'String', '  ... ', 'callback', @nope); %added by Jan 2018

switch upper(StimTypes{irec}),
    % Try to maintain as much uniformity amongst stimulus types as
    % possible.
    case {'RC','FS','TCKL','WAV','CAP'}
        set(h1, 'String', '1 PSTH', 'callback', @local_PSTH);
        set(h2, 'String', '2 cyclehisto', 'callback', @local_cyclehisto);
        set(h3, 'String', '3 ISI histo', 'callback', @local_FOISI);
        set(h4, 'String', '4 coeffvar', 'callback', @local_coeffvar)
        set(h5, 'String', '5 rateplot', 'callback', @local_rateplot);
        set(h6, 'String', '6 magn/phi', 'callback', @local_magn_phi); %added by Jan 2018

    case 'RF',
        set(h1, 'String', '1 PSTH', 'callback', @local_PSTH);
        set(h2, 'String', '2 IOfun', 'callback', @local_IOfun);
        set(h3, 'String', '3 area', 'callback', @local_resparea);
        set(h4, 'String', '4 coeffvar', 'callback', @local_coeffvar);
        set(h5, 'String', '5 rateplot', 'callback', @local_rateplot);
        set(h6, 'String', '6 average', 'callback', @local_anamean);
    
    case {'RCN','NITD','ARMIN','NRHO','MOVN','IRN','HP'}
        set(h1, 'String', '1 PSTH', 'callback', @local_PSTH);
        set(h2, 'String', '2 revcor', 'callback', @local_revcor);
        set(h3, 'String', '3 ISI histo', 'callback', @local_FOISI);
        set(h4, 'String', '4 magn/phi', 'callback', @local_magn_phi);
        set(h5, 'String', '5 rateplot', 'callback', @local_rateplot);
        
    case {'MTF','DEP','RCM','RAM','BBFC','BBFM','BBFB','ITD','ILD','MBL',...
            'CFS','CTD','CSPL','MOVING_W','ZWICKER','QFM','NSAM','HAR',...
            'ENH_W','ENH_D','ENH_DB','W','NPHI','ENH_2T','ENH_FC','ENH_DURC'},
        set(h1, 'String', '1 PSTH', 'callback', @local_PSTH);
        set(h2, 'String', '2 cyclehisto', 'callback', @local_cyclehisto);
        set(h3, 'String', '3 ISI histo', 'callback', @local_FOISI);
        set(h4, 'String', '4 magn/phi', 'callback', @local_magn_phi);
        set(h5, 'String', '5 rateplot', 'callback', @local_rateplot);
        set(h6, 'String', '6 average', 'callback', @local_anamean); %added by Jan 2018
        
    case 'MASK',
        set(h1, 'String', '1 rateplot', 'callback', @local_rateplot);
        set(h2, 'String', '2 maskcorr', 'callback', @local_maskcorr);
        
    case 'SUP',
        set(h1, 'String', '1 supspec', 'callback', @local_supspec);
        
    case 'ZW',
        set(h1, 'String', '1 apple', 'callback', @local_apple);
        
    case 'BINZW',
        set(h4, 'String', '4 rmsplot', 'callback', @local_rmsplot);
        
    case {'THR','CAP'}
        set(h1, 'String', '1 threshold curve', 'callback', @local_threshold_curve);
end % switch/case
%drawnow;

function local_custom_1(Src, Ev);
[figh, bt] = local_handle(Src, 'CustomAna_1_Button');
if ~isa(figh,'double')
    figh = figh.Number;
end;
feval(get(bt,'Callback'), Src, Ev);

function local_custom_2(Src, Ev);
[figh, bt] = local_handle(Src, 'CustomAna_2_Button');
if ~isa(figh,'double')
    figh = figh.Number;
end;
feval(get(bt,'Callback'), Src, Ev);

function local_custom_3(Src, Ev);
[figh, bt] = local_handle(Src, 'CustomAna_3_Button');
if ~isa(figh,'double')
    figh = figh.Number;
end;
feval(get(bt,'Callback'), Src, Ev);

function local_custom_4(Src, Ev);
[figh, bt] = local_handle(Src, 'CustomAna_4_Button');
if ~isa(figh,'double')
    figh = figh.Number;
end;
feval(get(bt,'Callback'), Src, Ev);

function local_custom_5(Src, Ev);
[figh, bt] = local_handle(Src, 'CustomAna_5_Button');
if ~isa(figh,'double')
    figh = figh.Number;
end;
feval(get(bt,'Callback'), Src, Ev);

function local_custom_6(Src, Ev);
[figh, bt] = local_handle(Src, 'CustomAna_6_Button');
if ~isa(figh,'double')
    figh = figh.Number;
end;
feval(get(bt,'Callback'), Src, Ev);

function local_keypress(Src, Ev);
%Ev
figh = parentfigh(Src);
switch Ev.Key,
    case 'escape';
        commandwindow;
    case 'delete';
        close(figh);
    case '1'; % 
        local_custom_1(figh, Ev);
    case '2'; % 
        local_custom_2(figh, Ev);
    case '3'; % 
        local_custom_3(figh, Ev);
    case '4'; % 
        local_custom_4(figh, Ev);
    case '5'; % 
        local_custom_5(figh, Ev);
    case '6'; % 
        local_custom_6(figh, Ev);
    case 'r'; % 
        local_dotraster(figh, Ev);
    case 'p'; % 
        local_paramview(figh, Ev);
    case 'd'; % 
        local_displaylist(figh, Ev);
    case 'u'; % 
        local_update(figh, Ev);
    case 't'; % 
        local_typelog(figh, Ev);
end
drawnow;
