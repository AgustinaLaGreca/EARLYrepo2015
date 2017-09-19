function ee=exist(E,N);
% experiment/exist - True for existing Experiment
%   exist(experiment(), 'Foo') - returns True if an Experiment named Foo 
%   exists on this computer, False otherwise. Existence is equivalent to the 
%   presence of the named subfolder in one of the Experimenter subfolders 
%   of parentfolder(experiment()). Note that the first input arg to exist
%   is a dummy experiment (see Methodizing).
%
%   See also Experiment/locate, Experiment/find.

if ~isvarname(N), % N must be varname to count as a valid experiment name
    ee = false;
else, % check existence of Foo.ExpDef in appropriate folder
    Fld = folder(dummy(experiment(), N));
    ee = ~isempty(Fld) && exist(fullfile(Fld, [N '.ExpDef']));
end








