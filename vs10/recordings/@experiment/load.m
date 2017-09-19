function E = load(E, Er);
% Experiment/load - load experiment from disk
%    E = load(E) re-loads E from disk to get the most recent version.
%
%    E = load(E, 'Pietje') launches UIgetdir to select an experiment done
%    by experimenter Pietje.
%
%    See also Experiment/find.

if nargin<2,
    warning('off', 'MATLAB:fieldConstructorMismatch');
    E = load([filename(E) '.ExpDef'], '-mat');
    warning('on', 'MATLAB:fieldConstructorMismatch');
    E = experiment(E.E);
else, % experimenter given: prompt user
    if ~isknownexperimenter(experiment, Er),
        error(['Experimenter ''' Er ''' not found on this computer.']);
    end
    Edir = fullfile(parentfolder(E), Er)
    Sf = subfolders(Edir)
    if isempty(Sf), error(['No experiments by ''' Er ''' found on this computer.']); end
    DD = fullfile(Edir, Sf{1})
    Nam = uigetdir(DD, 'Select Experiment');
    if isequal(0,Nam), % user cancelled
        E = experiment(); % void
    else,
        [dum, Nam] = fileparts(Nam);
        E = find(experiment, Nam);
    end
end









