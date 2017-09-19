function S = listall(dum, Experimenter); %#ok<INUSL>
% Experiment/listall - list all experiments stored on this computer
%    listall(experiment()) lists all experiments stored on this computer.
%    Note that the first argument is a dummy for Methodizing this function.
%    Listall does a recursive search in the subfolders of
%    parentfolder(experiment). 
%
%    listall(experiment,'Pietje') restricts the list to experiments by
%    experimenter Pietje.
%
%    S = experimenter(...) suppresses listing to the screen and returns a
%    struct S with fields Name, Experimenter, Type, Started, Modified,
%    StatusModified, Computer, and Location.
%
%    See Datadir, Experiment, DSread.

if nargin<2, 
    Ee = allexperimenters(dum);
else,
    if ~ismember(lower(Experimenter), lower(allexperimenters(dum))),
        warning(['Experimenter ''' Experimenter ''' not known on this computer.']);
    end
    Ee = cellify(Experimenter);
end

iexp=0;
for ii=1:numel(Ee),
    ErDir = fullfile(parentfolder(dum), Ee{ii});
    XP = subfolders(ErDir); % experiments of experimenter Ee{ii}
    for jj=1:numel(XP),
        iexp = iexp + 1;
        Exp = find(dum, XP{jj});
        S(iexp) = Exp.ID; % the ID field has it all
    end
end

if nargout<1,
    if exist('S', 'var'),
        disp(['Experiments stored on ' CompuName ':']);
        prettyprint({S.Name});
        clear S;
    else,
        disp(['No Experiments found.']);
    end
    disp(' ');
elseif ~exist('S', 'var'), % no experiments found. Still return a struct having the right fields.
    S = dum.ID;
    S = S([]);
end





