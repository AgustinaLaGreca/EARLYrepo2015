function ds = add_comment( ds, comment )
%add_comment: adds a comment to the ID struct of the dataset
%   ds: a dataset object
%   comment: a string containing the comment which will be added to de
%            comment field of the ID struct

if isempty(ds)
   error('add_comment.m :: the dataset is empty. Please supply a correct dataset'); 
end

if ischar(comment)
    ds.ID.comment = comment;
else
    warning('add_comment.m :: comment not added to the datast since it is not a string');
end


end

