function bh = paramview(D);
% dataset/paramview - view the stimulus parameters of a dataset
%    paramview(D) displays the stimulus parameters of dataset D in a GUI
%    identical to the one used at recording time.
%
%    See also stimGUI.

if numel(D)>1,
    for id=1:numel(D),
        paramview(D(id));
        set(gcg, 'position', get(gcg,'position')+[20*id, -20*id, 0, 0]);
    end
    return;
end

% single D from here

EXP = D.Stim.Experiment;
StimName = D.Stim.StimType;
Params = feval(fhandle(['stimdef' StimName]), EXP); % Params = stimdefXXX(EXP);

SG=GUIpiece([D.Stim.StimType 'Menu'],[],[0 0],[10 4]);
SG = add(SG, Params);
SG=marginalize(SG,[0 20]);

WinTitle = ['Stimulus parameters of ' name(EXP) ', rec ' num2str(D.ID.iDataset)  ' <' IDstring(D) '>'];
figh = newGUI([StimName 'menu'], WinTitle, {fhandle(mfilename), StimName, EXP});
%setGUIdata(figh, 'StimulusType', StimName);
draw(figh, SG); 
GUIfill(figh, D.Stim.GUIgrab);
% show play time, including baseline (if any)
BL = getFieldOrDefault(D.Stim.Presentation,'baseline' ,[]);
if isempty(BL), BLdur = 0;
else, BLdur = BL.PreDur + BL.PostDur;
end
Pres = D.Stim.Presentation;
meanISI  = mean(Pres.PresDur(2:end-1)); % ignore baselines
ReportPlayTime(figh, Pres.Ncond, Pres.Nrep, meanISI, BLdur);
% set enabled edit fields to inactive 
Q = getGUIdata(figh, 'Query');
eh = edithandle(Q);
en = get(eh,'enable');
ien = strmatch('on', en);
set(eh(ien),'enable', 'inactive');
% disable all buttons
bh = findobj(gcg,'type', 'uicontrol', 'style', 'pushbutton');
set(bh, 'enable', 'inactive', 'callback', [], 'ButtonDownFcn' ,[]);



