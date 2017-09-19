function F = filename(E);
% experiment/filename - full file name of experiment data on this computer
%    filename(E) returns the full filename prefix of teh experimental data
%    on the current computer, which is 
%        <datadir>/<experimenter>/<expname>/<expname> 
%    with no extension.
%    A list of all data files belonging to E is produced by
%    dir([filename(E) '*.*']). Per definition, Filename returns '' for a
%    void E.
%
%    See also Datadir, experiment/filename, experiment/save.

if isvoid(E), 
    F = ''; 
else,
    F = fullfile(folder(E), name(E));
end












