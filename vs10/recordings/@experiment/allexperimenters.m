function L = allexperimenters(E, Xnam);
% experiment/allexperimenters - list of all known experimenters
%   allexperimenters(experiment()) returns a list (cell string) of all known 
%   experimenters at this computer. This is simply a list of the names of
%   all subdirectories of parentfolder(experiment()).
%
%   Note that the first input arg is a dummy for Methodizing this function.
%   
%   See also Experiment/allexperimenters, Experiment/checkin, Methodizing.

L = subfolders(parentfolder(E));
for ii=1:numel(L),
    L{ii}(1) = upper(L{ii}(1)); % proper names
end






