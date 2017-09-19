function H = listall(dum);
% listall - list all hoard objects
%    listall(hoard) lists all hoard objects.
%
%    See also hoard, hoard/clearall, hoard/find, hoard/findall.

clear(dum); % dum serves to make this a hoard method. Thanks & bye.

hmarks = findall(0,'tag','I_am_a_lonesome_hoard');
for ii=1:numel(hmarks),
    hfig = get(hmarks(ii), 'parent');
    disp(hoard(hfig));
end








