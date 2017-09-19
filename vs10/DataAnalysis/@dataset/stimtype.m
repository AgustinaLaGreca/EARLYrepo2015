function T=stimtype(D);
% Dataset/stimtype - stimulus type
%   stimtype(D) returns a char string describing the stimulus type of D.
%
%   For arrays, a cell string is returned having the same size as D.
%
%   See also dataset/IDstring.

if numel(D)==1,
    T = upper(D.Stim.StimType);
else,
    sD = size(D);
    D = D(:).';
    T = [{stimtype(D(1))}, stimtype(D(2:end))];
    T = reshape(T, sD);
end



