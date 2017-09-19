function H = anaview(DS, Chan, iCond);
% Dataset/anaview - waterfall plot showing raw analog data
%    anaview(DS, Chan) plots the raw recordings in analog channel Chan
%    of dataset DS. (See dataset/anachan for channel specification.)
%    The individual reps are superimposed and the mean is plotted on top of
%    them usinga thicker line. Anaview returns the handles to the figure(s)
%    it opens.
%
%    D = anaview(DS, Chan, iCond) only selects the stimulus conditions
%    with indices in array iCond. 
%
%    See also Dataset/recview, Dataset/anadata, Dataset/anachan, 
%    Dataset/anamean.

MAX_NSUBPLOT = 10; % max # subplots in one fig

[iCond] = arginDefaults('iCond', 0); % default: all conditions
iCond = replaceMatch(iCond, {0 ':'}, 1:Ncond(DS));

NC = numel(iCond);
Nsubplot = min(NC, MAX_NSUBPLOT);
[ScaleFac, Unit] = conversionfactor(anachan(DS, Chan));
StimPres = DS.Stim.Presentation;

ifig = 0; 
H = [];
isubplot = Nsubplot; % force opening new fig
DfigPos = [0.017 -0.022 0 0];
Yoffset = 0;
[Ymin, Ymax] = deal(inf, 0);
for icond=1:NC,
    isubplot = isubplot + 1;
    if isubplot>Nsubplot, % open new figure, start with first subplot
        if ifig>0, % finish figure
            set(findobj(figh), 'visib', 'on');
            for ii=1:numel(TXT_left),
                text(min(xlim)+0.002*diff(xlim), TXT_left(ii).Y, TXT_left(ii).Str, 'color', 'k', 'fontweight', 'bold');
                text(max(xlim)+0.02*diff(xlim), TXT_right(ii).Y, TXT_right(ii).Str, 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
            end
            ylim([Ymin Ymax]);
            clear TXT_left TXT_right
        end
        figh = figure;
        ifig = ifig+1;
        H(ifig) = figh;
        isubplot = 1;
        [Ymin, Ymax] = deal(inf, 0);
        Yoffset = 0;
    end
    set(figh,'units', 'normalized', 'position', [0.165 0.188 0.684 0.727]+(numel(H)-1)*DfigPos);
    [D, dt, t0, DataType] = anadata(DS, Chan, iCond(icond));
    D = D*ScaleFac;
    meanD = mean(D(:));
    D = D-meanD;
    Yspacing = 5*mean(std(D,[],2));
    Yoffset = Yoffset - Yspacing;
    Doffset{ifig}(isubplot) = Yoffset;
    xdplot(dt, D+Yoffset, 'visib', 'off');
    xdplot(dt, mean(D,2)+Yoffset, 'w', 'linewidth', 2.5, 'visib', 'off');
    xdplot(dt, mean(D,2)+Yoffset, 'k', 'linewidth', 1.5, 'visib', 'off');
    TXT_left(isubplot) = struct('Y', Yoffset+0.2*Yspacing, 'Str', num2str(deciRound(meanD,2)));
    TXT_right(isubplot) = struct('Y', Yoffset, 'Str', [num2str(iCond(icond)), ': ' CondLabel(StimPres, iCond(icond))]);
%     text(max(xlim)+0.02*diff(xlim), Yoffset, [num2str(iCond(icond)), ': ' CondLabel(StimPres, icond)]);
%     text(min(xlim)+0.002*diff(xlim), Yoffset+0.2*Yspacing, num2str(deciround(meanD,2)), 'fontsize', 8, 'fontweight', 'bold', 'color', 0.05*[1 1 1]);
    Yoffset = Yoffset - Yspacing;
    Ymin = min(Ymin, Yoffset);
end
% finish last fig
set(findobj(figh), 'visib', 'on');
for ii=1:numel(TXT_left),
    text(min(xlim)+0.002*diff(xlim), TXT_left(ii).Y, TXT_left(ii).Str, 'color', 'k', 'fontweight', 'bold');
    text(max(xlim)+0.02*diff(xlim), TXT_right(ii).Y, TXT_right(ii).Str, 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
end
if NC>1,
    ylim([Ymin Ymax]);
end






