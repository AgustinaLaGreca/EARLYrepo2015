function cleanup(G, varargout);
% comport_data/cleanup - cleanup temporary data associated with comport_data.
%   cleanup(D) removes any temporary data associated with D from disk.
%   To date, comport_data objects do not have temporary data on disk, so
%   this method is void. It is needed anyhow, because cleanup must be
%   overloaded for all grabbeddata subclasses.
%
%   See also grabbedata/cleanup, dataset/cleanup.

eval(IamAt);







