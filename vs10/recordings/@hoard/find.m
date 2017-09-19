function H = find(dum, Name);
% find - find named hoard object
%    find(hoard, 'Foo') returns hoard object named Foo, if any.
%
%    See also hoard, hoard/findall, hoard/listall.

clear(dum); % dum serves to make this a hoard method. Thanks & bye.

if isnumeric(Name), Name = num2str(Name); end

hmarks = findall(0,'tag','I_am_a_lonesome_hoard');
H = [];
for ii=1:numel(hmarks),
    hfig = get(hmarks(ii), 'parent');
    if isequal(Name, get(hfig,'tag')),
        H = [H hoard(hfig)];
    end
end








