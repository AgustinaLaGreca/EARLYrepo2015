function StimPlot(S);
% StimPlot - plot stimulus waveforms
%   StimPlot(S), where S is the stimulus-defining struct (as in the Stim
%   field of dataset objects) opens a GUI for viewing the stimulus
%   waveforms.
%   
%   StimPlot(h), where h is the graphics handle of a stimulus menu created
%   by StimGUI, uses the stimulus of that menu. 

% get stimuli from stimulus GUI if needed

S = local_getstim(S);
if isempty(S), return; end;
%-----------the whole GUI-----------------
G=GUIpiece('StimPlot',[],[0 0],[10 4]);
if strcmpi(S.StimType,'hp') || strcmpi(S.StimType,'armin')|| (strcmpi(S.StimType,'har') && length(S.Waveform(1,:)) ~= 1) || strcmpi(S.StimType,'irn') || strcmpi(S.StimType,'qfm')
    AXWL = AxesDisplay('@WaveformAxesL', [580 280], 'linewidth', 2, 'fontsize',12); % Left
    AXML = AxesDisplay('MagnitudeAxesL', [580 280], 'linewidth', 2, 'fontsize',12);
    AXPL = AxesDisplay('PhaseAxesL', [580 280], 'linewidth', 2, 'fontsize',12);
    AXWD = AxesDisplay('@WaveformAxesD', [580 280], 'linewidth', 2, 'fontsize',12); % Difference
    AXMD = AxesDisplay('MagnitudeAxesD', [580 280], 'linewidth', 2, 'fontsize',12);
    AXPD = AxesDisplay('PhaseAxesD', [580 280], 'linewidth', 2, 'fontsize',12);
    AXWR = AxesDisplay('@WaveformAxesR', [580 280], 'linewidth', 2, 'fontsize',12); % Right
    AXMR = AxesDisplay('MagnitudeAxesR', [580 280], 'linewidth', 2, 'fontsize',12);
    AXPR = AxesDisplay('PhaseAxesR', [580 280], 'linewidth', 2, 'fontsize',12);
    UP = local_userPanel(S);
    G = add(G, UP);
    G = add(G, AXWL, nextto(UP), [-200 0]);
    G = add(G, AXML, nextto(UP), [-200 220]);
    G = add(G, AXPL, nextto(UP), [-200 440]);
    G = add(G, AXWD, nextto(UP), [125 0]);
    G = add(G, AXMD, nextto(UP), [125 220]);
    G = add(G, AXPD, nextto(UP), [125 440]);
    G = add(G, AXWR, nextto(UP), [475 0]);
    G = add(G, AXMR, nextto(UP), [475 220]);
    G = add(G, AXPR, nextto(UP), [475 440]);
else
    AXW = AxesDisplay('@WaveformAxes', [900 280], 'linewidth', 2, 'fontsize',12);
    AXM = AxesDisplay('MagnitudeAxes', [900 280], 'linewidth', 2, 'fontsize',12);
    AXP = AxesDisplay('PhaseAxes', [900 280], 'linewidth', 2, 'fontsize',12);
    UP = local_userPanel(S);
    G = add(G, UP);
    G = add(G, AXW, nextto(UP), [-100 0]);
    G = add(G, AXM, nextto(UP), [-100 220]);
    G = add(G, AXP, nextto(UP), [-100 440]);    
end
% SG = add(SG, Pedit);
% SG = add(SG, Pview);
% SG = add(SG, Params);
% SG = add(SG, local_Action, 'below', [20 10]);
% SG = add(SG, local_Dataview, 'nextto', [10 0]);
G=marginalize(G,[0 0]);
% draw the GUI
figh = newGUI('Stimplot', [S.StimType ' stimulus plot'], {fhandle(mfilename), S});
setGUIdata(figh, 'Stim', S);
draw(figh, G); 
TracePlotInterface(figh,false); % don't link axes
set(figh, 'toolbar','figure');
Q = getGUIdata(figh,'Query');
set(edithandle(Q('iCond')),'string', '1');
local_refresh(figh); % realize first plot

%=====================
function S = local_getstim(S);
% get stim
if isGUI(S),
    if ~StimCheck(S), 
        S = [];
        return; 
    end;
    S = getGUIdata(S, 'StimParam', []);
end

function UP = local_userPanel(S);
UP = GUIpanel('Select', 'Select');
Pres = S.Presentation;
[X, Y] = deal(Pres.X, Pres.Y);
Ncond = Pres.Ncond;
iCond = ParamQuery('iCond', 'Stim #', 'XXXXX', '', ...
    'rreal/posint', ['Index of stimulus condition. 1...' num2str(Ncond) '. Hit <return> to activate.'],1, 'callback', @local_refresh);
StimInfoX=messenger('StimInfo_X', sprintf([X.ParName ': ' X.FormatString], max(X.PlotVal)),1, ...
    'fontweight', 'bold', 'ForegroundColor', [0 0 0.7]);
if ~isempty(Y),
    StimInfoY=messenger('StimInfo_Y', sprintf([Y.ParName ': ' Y.FormatString], max(Y.PlotVal)),1, ...
        'fontweight', 'bold', 'ForegroundColor', [0.7 0 0]);
end
MM=messenger('@MainMess', 'XXXXXXXXXXXXXXXXXXXXXXXXXX', 4, ...
    'fontweight', 'bold', 'fontsize',13, 'backgroundcolor', [0.8 1 0.8]);
RF = ActionButton('Recompute', 'Recompute', 'Recompute', ...
    'Recompute stimulus using the current stimulus parameters from the stimulus GUI.',...
    @local_recompute);
HLP = ActionButton('Help', 'Help', 'XXXXXXX', ...
    'Click to display help on fast keys for waveform plot navigation.',...
    @local_help);
UP = add(UP, iCond);
UP = add(UP, StimInfoX, 'below');
if ~isempty(Y),
    UP = add(UP, StimInfoY, 'below');
end
UP = add(UP, RF, 'below');
UP = add(UP, HLP, nextto(RF), [80 0]);
UP = add(UP, MM, below(RF), [0 20]);
UP = marginalize(UP, [10 250]);

function local_refresh(Src, varargin);
% refresh plot
figh = parentfigh(Src);
figh = figh.Number;
GV = GUIval(figh);
if isempty(GV), return; end
% get stimulus from GUI & plot it
S = getGUIdata(figh, 'Stim');
Pres = S.Presentation;
if GV.iCond>Pres.Ncond,
    GUImessage(figh,['Index exceeds # conditions = ' num2str(Pres.Ncond) '.'], 'error', 'iCond');
    return;
end
GUImessage(figh,' ');
W = S.Waveform(GV.iCond, :);
[X, Y] = deal(Pres.X, Pres.Y);
set(0,'showhidden', 'on');
% plot waveform; for pitch stimuli plot for each ear and plot difference
if strcmpi(S.StimType,'hp') || strcmpi(S.StimType,'armin') || (strcmpi(S.StimType,'har') && length(S.Waveform(1,:)) ~= 1) || strcmpi(S.StimType,'irn') || strcmpi(S.StimType,'qfm')
    ah = GUIaxes(figh, '@WaveformAxesL');
    axes(ah);
    if strcmpi(W(1).DAchan,'l')
        channel_index = 1;
    else
        channel_index = 2;
    end
    plot(W(channel_index));
    Atten = S.Attenuation.AnaAtten;
    AttStr = ['Analog attenuation: (' trimspace(num2str(0.1*round(10*Atten))), ') dB'];
    title(AttStr,'FontWeight','bold');
    set(gca,'ylimmode', 'auto');
    ylim(ylim); % fix y limits so they won't jump when scrolling
    set(gca,'ylimmode', 'manual');
    % plot magnitude spectrum
    ah = GUIaxes(figh, 'MagnitudeAxesL');
    axes(ah);
    plotMagnitude(W(channel_index));
    title('Magnitude','FontWeight','bold');
    set(gca,'ylimmode', 'auto');
    ylim(ylim); % fix y limits so they won't jump when scrolling
    if isfield(S,'HighFreq')
        xlim([0 S.HighFreq*1.1]);
    elseif isfield(S,'Fc')
        xlim([0 S.Fc(end)*1.1]);
    else
        xlim([0 max(S.Fcar*max(max([S.F01Harmonics, S.F02F01.*S.F02Harmonics])))*1.1]); 
    end
    set(gca,'ylimmode', 'manual');
    % plot phase spectrum
    ah = GUIaxes(figh, 'PhaseAxesL');
    axes(ah);
    [ph,fr] = plotPhase(W(channel_index));
    title('Phase','FontWeight','bold');
    set(gca,'ylimmode', 'auto'); 
    ylim(ylim); % fix y limits so they won't jump when scrolling
    if isfield(S,'HighFreq')
        xlim([0 S.HighFreq*1.1]);
    else
       xlim([0 max(S.Fcar*max(max([S.F01Harmonics, S.F02F01.*S.F02Harmonics])))*1.1]); 
    end
    set(gca,'ylimmode', 'manual');
    %set(0,'showhidden', 'off');
    Xm = find(messenger,figh, 'StimInfo_X');
    Xmess = sprintf([X.ParName ': ' X.FormatString], X.PlotVal(GV.iCond));
    report(Xm, Xmess);
    if ~isempty(Y),
        Ym = find(messenger,figh, 'StimInfo_Y');
        Ymess = sprintf([Y.ParName ': ' Y.FormatString], Y.PlotVal(GV.iCond));
        report(Ym, Ymess);
    end
    
    %Right
    
    ah = GUIaxes(figh, '@WaveformAxesR');
    axes(ah);
    if strcmpi(W(1).DAchan,'r')
        channel_index = 1;
    else
        channel_index = 2;
    end
    plot(W(channel_index));
    Atten = S.Attenuation.AnaAtten;
    AttStr = ['Analog attenuation: (' trimspace(num2str(0.1*round(10*Atten))), ') dB'];
    title(AttStr,'FontWeight','bold');
    set(gca,'ylimmode', 'auto');
    ylim(ylim); % fix y limits so they won't jump when scrolling
    set(gca,'ylimmode', 'manual');
    % plot magnitude spectrum
    ah = GUIaxes(figh, 'MagnitudeAxesR');
    axes(ah);
    plotMagnitude(W(channel_index));
    title('Magnitude','FontWeight','bold');
    set(gca,'ylimmode', 'auto');
    ylim(ylim); % fix y limits so they won't jump when scrolling
    if isfield(S,'HighFreq')
        xlim([0 S.HighFreq*1.1]);
    else
       xlim([0 max(S.Fcar*max(max([S.F01Harmonics, S.F02F01.*S.F02Harmonics])))*1.1]); 
    end
    set(gca,'ylimmode', 'manual');
    % plot phase spectrum
    ah = GUIaxes(figh, 'PhaseAxesR');
    axes(ah);
    [ph,fr] = plotPhase(W(channel_index));
    title('Phase','FontWeight','bold');
    set(gca,'ylimmode', 'auto'); 
    ylim(ylim);% fix y limits so they won't jump when scrolling
    if isfield(S,'HighFreq')
        xlim([0 S.HighFreq*1.1]);
    else
       xlim([0 max(S.Fcar*max(max([S.F01Harmonics, S.F02F01.*S.F02Harmonics])))*1.1]); 
    end
    set(gca,'ylimmode', 'manual');
    %set(0,'showhidden', 'off');
    Xm = find(messenger,figh, 'StimInfo_X');
    Xmess = sprintf([X.ParName ': ' X.FormatString], X.PlotVal(GV.iCond));
    report(Xm, Xmess);
    if ~isempty(Y),
        Ym = find(messenger,figh, 'StimInfo_Y');
        Ymess = sprintf([Y.ParName ': ' Y.FormatString], Y.PlotVal(GV.iCond));
        report(Ym, Ymess);
    end
    
    % Difference
    ah = GUIaxes(figh, '@WaveformAxesD');
    axes(ah);
    if strcmpi(W(1).DAchan,'r')
        channel_index = 1;
    else
        channel_index = 2;
    end
    plot(W,'diff');
    Atten = S.Attenuation.AnaAtten;
    AttStr = ['Analog attenuation: (' trimspace(num2str(0.1*round(10*Atten))), ') dB'];
    title(AttStr,'FontWeight','bold');
    set(gca,'ylimmode', 'auto');
    ylim(ylim); % fix y limits so they won't jump when scrolling
    set(gca,'ylimmode', 'manual');
    % plot magnitude spectrum
    ah = GUIaxes(figh, 'MagnitudeAxesD');
    axes(ah);
    plotMagnitude(W,'diff');
    title('Magnitude','FontWeight','bold');
    set(gca,'ylimmode', 'auto');
    ylim(ylim); % fix y limits so they won't jump when scrolling
    if isfield(S,'HighFreq')
        xlim([0 S.HighFreq*1.1]);
    else
       xlim([0 max(S.Fcar*max(max([S.F01Harmonics, S.F02F01.*S.F02Harmonics])))*1.1]); 
    end
    set(gca,'ylimmode', 'manual');
    % plot phase spectrum
    ah = GUIaxes(figh, 'PhaseAxesD');
    axes(ah);
    [ph,fr] = plotPhase(W,'diff');
    title('Phase Diff (R-L)','FontWeight','bold');
    set(gca,'ylimmode', 'auto'); 
    if isfield(S,'HighFreq')
        if ph(round(S.HighFreq*1.1/fr)) < 0
            ylim([ph(round(S.HighFreq*1.1/fr)) 0]); % fix y limits so they won't jump when scrolling
        elseif ph(round(S.HighFreq*1.1/fr)) > 0
            ylim([min(ph(1:round(S.HighFreq*1.1/fr))) max(ph(1:round(S.HighFreq*1.1/fr)))]);
        end
    else
        fre=max(S.Fcar*max(max([S.F01Harmonics, S.F02F01.*S.F02Harmonics])));
        if ph(round(fre*1.1/fr)) < 0
            ylim([ph(round(fre*1.1/fr)) 0]); % fix y limits so they won't jump when scrolling
        elseif ph(round(fre*1.1/fr)) > 0
            ylim([min(ph(1:round(fre*1.1/fr))) max(ph(1:round(fre*1.1/fr)))]);
        end
        
    end
        
    if isfield(S,'HighFreq')
        xlim([0 S.HighFreq*1.1]);
    else
       xlim([0 max(S.Fcar*max(max([S.F01Harmonics, S.F02F01.*S.F02Harmonics])))*1.1]); 
    end
    set(gca,'ylimmode', 'manual');
    %set(0,'showhidden', 'off');
    Xm = find(messenger,figh, 'StimInfo_X');
    Xmess = sprintf([X.ParName ': ' X.FormatString], X.PlotVal(GV.iCond));
    report(Xm, Xmess);
    if ~isempty(Y),
        Ym = find(messenger,figh, 'StimInfo_Y');
        Ymess = sprintf([Y.ParName ': ' Y.FormatString], Y.PlotVal(GV.iCond));
        report(Ym, Ymess);
    end
else
    ah = GUIaxes(figh, '@WaveformAxes');
    axes(ah);
    plot(W);
    Atten = S.Attenuation.AnaAtten;
    AttStr = ['Analog attenuation: (' trimspace(num2str(0.1*round(10*Atten))), ') dB'];
    title(AttStr,'FontWeight','bold');
    set(gca,'ylimmode', 'auto');
    ylim(ylim); % fix y limits so they won't jump when scrolling
    set(gca,'ylimmode', 'manual');
    % plot magnitude spectrum
    ah = GUIaxes(figh, 'MagnitudeAxes');
    axes(ah);
    plotMagnitude(W);
    title('Magnitude','FontWeight','bold');
    set(gca,'ylimmode', 'auto');
    ylim(ylim); % fix y limits so they won't jump when scrolling
%     xlim([0 S.Fcar(end)*1.1]);
    set(gca,'ylimmode', 'manual');
    % plot phase spectrum
    ah = GUIaxes(figh, 'PhaseAxes');
    axes(ah);
    plotPhase(W);
    title('Phase','FontWeight','bold');
    set(gca,'ylimmode', 'auto'); 
    ylim(ylim); % fix y limits so they won't jump when scrolling
%     xlim([0 S.Fcar(end)*1.1]);
    set(gca,'ylimmode', 'manual');
    %set(0,'showhidden', 'off');
    Xm = find(messenger,figh, 'StimInfo_X');
    Xmess = sprintf([X.ParName ': ' X.FormatString], X.PlotVal(GV.iCond));
    report(Xm, Xmess);
    if ~isempty(Y),
        Ym = find(messenger,figh, 'StimInfo_Y');
        Ymess = sprintf([Y.ParName ': ' Y.FormatString], Y.PlotVal(GV.iCond));
        report(Ym, Ymess);
    end
end
function local_recompute(Src,Ev,LR);
figh = parentfigh(Src); % stimplot figure
figh = figh.Number;
S = getGUIdata(figh, 'Stim');
hstim = S.handle.GUIfig; % stimulus GUI
if ~isGUI(S.handle.GUIfig),
    GUImessage(figh, 'Cannot find stimulus GUI.', 'error');
    return;
end
if ~StimCheck(hstim), return; end
S = getGUIdata(hstim, 'StimParam');
setGUIdata(figh, 'Stim', S);
local_refresh(figh);
% bug: remove redundant figure
if ishandle(1), close(1); end

function local_help(Src, Ev, LR);
GUImessage(parentfigh(Src), 'See command window');
help('TracePlotInterface');
commandwindow;



