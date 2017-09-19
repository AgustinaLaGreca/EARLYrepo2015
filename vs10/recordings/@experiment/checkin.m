function E=checkin(E, Nam, Er, T);
% experiment/checkin - register Experiment object
%   E=checkin(E, 'RG101999', 'Bas', 'ZWOAE') registers experiment E having
%   name RG101999, experiment Bas and experiment type ZWOAE. Its directory
%   is created, the (yet unfilled) E is saved, and E is made the current
%   experiment. The name must be unique on this computer. The Experimenter
%   must be known (that is, his/her subfolder of parentfolder()).

ID.Name = Nam;
ID.Experimenter = username();
ID.Type = T;

% name must be unique (across all experimenters!)
if ~isvarname(Nam), 
    error('Experiment name must be valid Matlab identifier (see ISVARNAME).');
elseif exist(experiment, Nam), 
    error(['Experiment name ''' Nam ''' has been used before (' locate(experiment, Nam, 1) '). Use a unique name.']);
elseif ~isknownexperimenter(E, Er),
    error(['Experimenter ''' Er ''' not known on this computer. See Experiment/allexperimenters.']);
end

ID.Started = datestr(now);
ID.Computer = CompuName;
ID.Location = Where;

if ~isequal(fieldnames(E.ID), fieldnames(ID)),
    error('ID fieldnames of Experiment object are incompatible with experiment/checkin. Update checkin!.');
end
E.ID = ID;

% make this experiment the current one
makecurrent(E);

% Save
% create experiment directory
[okay, Mess] = mkdir(fullfile(parentfolder(E), 'Exp'), Nam);
error(Mess);
save(E,['Experiment ' name(E) ' defined ' E.ID.Started '.'], 'new');








