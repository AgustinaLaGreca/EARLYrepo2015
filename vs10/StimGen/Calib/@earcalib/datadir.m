function D = datadir(EC);
% Earcalib/datadir - get or set directory for Earcalib data
%   datadir(EC) returns the directory where the Earcalib object EC is saved. 
%   This is the data directory of the experiment to which EC belongs:
%   folder(experiment(EC)).
%   
%   See also Datadir, setuplist, Methodizing.

D = folder(EC.Experiment);





