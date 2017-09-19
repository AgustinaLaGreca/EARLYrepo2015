function T=stimtype(D);
% Pooled_dataset/stimtype - stimulus type
%   stimtype(D) returns a cell array of strings contining stimtype of its pooled members
%
%   For arrays, a cell string is returned having the same size as D.
%
%   See also dataset/IDstring.


T = cellify(stimtype(D.DS));



