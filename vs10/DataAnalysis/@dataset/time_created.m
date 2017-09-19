function C = time_created(DS);
% Dataset/time_created - creation time of dataset
%    time_created(DS) returns the time DS was created, i.e., DS.ID.created.
%    For arrays DS, a cell array is returned.
%    See also dataset.

if isa(DS, 'pooled_dataset'),
    DS = members(DS);
    ID = [DS.ID];
    for ii=1:numel(ID),
       dn(ii) = datenum((ID(ii).created));
    end
    [dum, imax] = max(dn);
    DS = DS(imax);
end
id = [DS(:).ID];
C = {id.created};
C = reshape(C, size(DS));
if numel(C)==1, C = C{1}; end








