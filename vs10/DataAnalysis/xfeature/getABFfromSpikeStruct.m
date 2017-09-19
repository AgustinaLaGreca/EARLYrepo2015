function [D, DS, EE] = getABFfromSpikeStruct(S);
% getABFfromSpikeStruct - retrieve ABF data that produced xspikes output S
%   usage: D = getABFfromSpikeStruct(S), where S is return arg of xspikes.
%
%  [D, DS, EE] = getABFfromSpikeStruct(S); also tries to get SGSR dataset
%  DS and entry E from TKpool database. This is only possible if the ori
%
%  See also xspikes, readABFdata, readTKABF.

if isempty(S.ExpID), % no SGSR info; attempt to retrieve from ABFfile directly
    D = readABFdata(S.ABFname);
    [DS, EE] = deal([]);
else,
    [D, DS, EE] = readTKABF(S.ExpID, S.RecID, S.icond);
end  







