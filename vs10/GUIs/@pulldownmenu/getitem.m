function [Item, imatch, Mess] = getitem(P, Name);
% PulldownMenu/finditem - find named item contained in PulldownMenu object
%   finditem(P, 'Foo') returns the item  named Foo that is contained in
%   PulldownMenu object P. If no such named item is found, [] is returned.
%   If multiple items have the same name Foo, [] is returned. 
%   For elementary items of P (see PulldownMenu/additem), the Label
%   property is used as name.
%
%   [Item, I, imatch] = getitem(P, 'Foo'); also returns a message Mess
%   explaining why [] was returned, or empty if a unique Foo was found,
%   and index I of item Foo within the items contained in P.

%   See also PulldownMenu, PulldownMenu/additem.

Item = [];
Mess = '';
% different items in P may have different classes (see
% PulldownMenu/additem) so use a for loop.
imatch = [];
for ii=1:numel(P.Items),
    it = P.Items{ii};
    if isstruct(it), it.Name=it.Label; end
    if isequal(Name, it.Name),
        imatch = [imatch ii];
    end
end

if isempty(imatch),
    Mess = ['No item named ''' Name ''' found in Pulldown Menu ''' P.Name '''.'];
elseif numel(imatch)>1,
    Mess = ['Multiple items named ''' Name ''' found in Pulldown Menu ''' P.Name '''.'];
else,
    Item = P.Items{ii};
end
