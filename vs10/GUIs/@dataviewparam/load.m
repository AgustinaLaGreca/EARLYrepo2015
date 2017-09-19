function P = load(Pdefault, DV, Filename, DefName);
% dataviewparam/load - load dataviewparam
%    load(Pdef, XXX, 'File') loads the viewdataparam object previously
%    saved by save(P, 'File'), where XXX is P's Dataviewer. If none is 
%    found, Pdef is returned.
%
%    load(Pdef, XXX, '?') prompts the user for the filename,
%    then loads the object unless the user cancelled.
%
%    load(Pdef, XXX, '?', 'Foo') prompts the user for the filename while
%    suggesting 'Foo' as a default name, then loads the object unless the 
%    user cancelled.
%
%    If the definition of the previously stored dataviewparam has become
%    obsolete due to changes in its dataviewer, any new parameters are
%    taken from the default dataviewparam having this Dataviewer; any
%    obsolete parameters are removed.
%
%    See also dataviewparam/save, dataviewparam/getdefault, dataset/dotraster.


if nargin<4, Subdir = {}; end

DV = char(DV);
sdir = savedir(Pdefault, DV);
if isequal('?', Filename),
    Filefilter = fullfile(sdir, '*.dvparam');
    if ~isempty(DefName), DefName = strrep(Filefilter, '*', DefName); end
    Filename = uigetfile(Filefilter, 'Load dataview parameters from file. ''Def'' is the default set.', DefName);
    if isequal(0, Filename), P=Pdefault; return; end
    Filename = strtok(Filename, '.'); % remove extension
end

FN = [Filename, '.dvparam'];
FFN = fullfile(sdir, FN);

P = getcache(FFN, 1);
if isequal([], P), % return first arg (see help)
    P = Pdefault;
    if ~isvoid(P) && ~isequal(char(P.Dataviewer), DV),
        error('Non-void Default dataviewerparam has the wrong Dataviewer.');
    end
else, % merge with newest default
    Pnew = getdefault(P);
    P.Param = union(Pnew.Param, P.Param);
    P.Param = structpart(P.Param, fieldnames(Pnew.Param));
end











