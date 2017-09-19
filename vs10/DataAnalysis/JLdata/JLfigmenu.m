function JLfigmenu(fh, S);
% JLfigmenu - create pulldown menu in figure for single-condition JLxxx analyses
%     JLfigmenu(fh, S) 

setGUIdata(fh, 'JLbeatOutput', S); 

hm = uimenu('label', '&AnaData');

CommonProps = {'Parent', hm, 'callback', @local_cb};
CommonPropsSep = {'Parent', hm, 'callback', @local_cb, 'Separator', 'on'};

uimenu('label', [S.TTT], 'Parent', hm);

uimenu('label', 'View Waveform', 'Tag', 'ViewWaveform', CommonPropsSep{:});
uimenu('label', 'Clampfit', 'Tag', 'Clampfit', CommonProps{:});
uimenu('label', 'JLspontAna', 'Tag', 'JLspontAna', CommonProps{:});
uimenu('label', 'JLspontPeak', 'Tag', 'JLspontPeak', CommonProps{:});

uimenu('label', 'Display data struct', 'Tag', 'Display', CommonPropsSep{:});

uimenu('label', 'JLbeatPlot (full)', 'Tag', 'JLbeatPlot', CommonPropsSep{:});
uimenu('label', 'JLbeatPlot (spec)', 'Tag', 'JLbeatPlot2', CommonPropsSep{:});
uimenu('label', 'JLbeatPlot (wav)', 'Tag', 'JLbeatPlot3', CommonProps{:});
uimenu('label', 'JLbeatPlot (lin)', 'Tag', 'JLbeatPlot4', CommonProps{:});
uimenu('label', 'JLbeatPlot (STA)', 'Tag', 'JLbeatPlot5', CommonProps{:});
uimenu('label', 'JLbeatPlot (CHist)', 'Tag', 'JLbeatPlot6', CommonProps{:});
uimenu('label', 'JLbeatPlot (APs)', 'Tag', 'JLbeatPlot7', CommonProps{:});
%
hmi = uimenu('label', 'JLbeatVar', 'Tag', 'JLbeatVar', CommonPropsSep{:});
hmi = uimenu('label', 'JLprincomp', 'Tag', 'JLprincomp', CommonProps{:});


% ==================
function local_cb(Src, Ev);
fh = gcbf;
S = getGUIdata(fh, 'JLbeatOutput');
switch get(Src, 'Tag'),
    case 'Clampfit', clampfit(S);
    case 'ViewWaveform', local_view(S);
    case 'JLspontAna', JLspontAna(S);
    case 'JLspontPeak', JLspontPeaks(S,[],[],1);

    case 'Display', disp(S);

    case 'JLbeatPlot', JLbeatPlot(S);
    case 'JLbeatPlot2', JLbeatPlot(S, 2);
    case 'JLbeatPlot3', JLbeatPlot(S, 3);
    case 'JLbeatPlot4', JLbeatPlot(S, 4);
    case 'JLbeatPlot5', JLbeatPlot(S, 5);
    case 'JLbeatPlot6', JLbeatPlot(S, 6);
    case 'JLbeatPlot7', JLbeatPlot(S, 7);
        
    case 'JLbeatVar', JLbeatVar(S);
    case 'JLprincomp', JLprincomp(S);
end

function local_view(S);
% view waveform
D = readTKABF(S);
Sam = D.AD(1,1).samples;
fh = figure; JLfigmenu(fh,S);
subplot(2,1,1);
set(gcf,'units', 'normalized', 'position', [0.00625 0.229 0.987 0.681]);
dplot(D.dt_ms, Sam);
fenceplot(S.SPTraw, ylim, 'r');
subplot(2,1,2);
Sam = smoothen(diff(Sam), 5);
dplot([D.dt_ms D.dt_ms/2], Sam);
fenceplot(S.SPTraw, ylim, 'r');
TracePlotInterface(fh);





