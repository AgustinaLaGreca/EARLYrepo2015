function cleanup(D);
% Dataset/cleanup - remove temporary data on disk associated with dataset
%    cleanup(DS) checks whether any temporary data ("hoards") are
%    associated with DS. If so, these are deleted. Not intended for direct
%    use; called by dataset/save when the recording of DS was interrupted
%    and the user does not want to save the data.
%
%    The real work is delegated to datagrabber/cleanup.
%
%    See Dataset/save.

eval(IamAt);

Dtypes = fieldnames(D.Data); 
for ii=1:numel(Dtypes), % visit the D.Data fields & get them data
    dtype = Dtypes{ii};
    datafield = D.Data.(dtype); % current value; may either hold Grabbeddata or Datagrabber
    if isa(datafield, 'address'), datafield = get(datafield); end % pointer -> content
    if isa(datafield,'grabbeddata'), % might need cleanup; delegate to supplier
        cleanup(datafield);
    end
end


