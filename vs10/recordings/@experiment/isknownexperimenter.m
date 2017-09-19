function iv=isknownexperimenter(E, Xnam);
% experiment/isknownexperimenter - True for known experimenter
%   isknownexperimenter(E,'Foo') returns True when Foo is a known 
%   experimenter, Flase otherwise. Experimenters are known by their having
%   a named subdirectory of the parentfolder(Experiment()) directory.
%
%   Note that the first input arg is a dummy for Methodizing this function.
%   
%   See also Experiment/allexperimenters, Experiment/checkin, Methodizing.

iv = ismember(lower(Xnam), lower(allexperimenters(E)));









