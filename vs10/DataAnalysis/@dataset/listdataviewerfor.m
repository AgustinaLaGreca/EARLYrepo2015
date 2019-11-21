function DVlist = listdataviewerfor(dum, StimType)
% Dataset/listdataviewerfor - list of dataviewers for dataset analysis
%    ListDataviewer(dataset()) returns the names of dataset methods that
%    are currently registered as "dataviewers" for stimulus type stim.
%       
%    See also dataset/listdataviewer, databrowse.

%CHANGELOG
%by Abel: 10 may 2016: Add rateplot to RCM list of viewers


DVlist = listdataviewer(dum, 'get'); % intersect with this so we won't make "hard-coded" mistakes

% For now, the "hard-coded" part is copy paste from databrowse...
switch upper(StimType)
    case {'RC','FS','TCKL','WAV'}
        DVlist = intersect(DVlist,...
            {'dotraster', 'PSTH', 'FOISI', 'AOISI', 'coeffvar','rateplot','cyclehisto'});
    case {'RF','RAM'}
        DVlist = intersect(DVlist,...
            {'dotraster', 'PSTH', 'FOISI', 'AOISI', 'coeffvar','resparea','rateplot', 'cyclehisto'});
    case {'RCN','NITD','ARMIN','NRHO','MOVN','IRN','HP','BPN'}
        DVlist = intersect(DVlist,...
            {'dotraster', 'PSTH', 'FOISI', 'AOISI', 'revcor','rateplot', 'cyclehisto','AutoCorr'});
    case {'MTF','DEP','BBFC','BBFM','BBFB','ITD','ILD','MBL',...
            'CFS','CSPL','MOVING_W','ZWICKER','QFM','NSAM',...
            'ENH_W','ENH_D','ENH_DB','W','NPHI','ENH_2T','ENH_FC','ENH_DURC',...
            'SUP','ZW','BINZW'},
        DVlist = intersect(DVlist,...
            {'dotraster', 'PSTH', 'FOISI', 'AOISI', 'cyclehisto','rateplot'});
    case {'MASK','RCM','CTD'}
        DVlist = intersect(DVlist,...
            {'dotraster', 'PSTH', 'FOISI', 'AOISI', 'rateplot','cyclehisto'});
    case {'THR','CAP'}
        DVlist = {};
    case {'HAR','HARHAR'}
        DVlist = intersect(DVlist,...
            {'dotraster', 'PSTH', 'FOISI', 'AOISI', 'cyclehisto','rateplot'});
end















