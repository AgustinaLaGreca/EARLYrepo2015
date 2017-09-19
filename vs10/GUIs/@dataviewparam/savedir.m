function d = savedir(P, DV);
% dataviewparam/savedir - path for saving dataviewparam object.
%   savedir(P) returns the full path of the directory where P will be
%   stored when calling save(P). This directory depends on the dataviewer
%   of P. The directory will be created if it does not exist yet.
%
%   save(Pvoid, DV), where Pdum is a void dataviewparam object, and DV is a
%   dataviewer, returns the directory where dataviewparam object haviung
%   dataviewer DV will be stored.
%
%   See also dataviewparam/save, dataviewparam/load, dataset/isdataviewer.

if nargin<2, 
    DV = dataviewer(P);
end
DV = lower(char(DV));
if ~isdataviewer(dataset(),DV),
    error(['''' DV ''' is not an existing dataviewer method. See dataset/isdataviewer.']);
end
d = GUIdefaultsDir(['DataviewParam\' DV], 'create'); % create if necessary





