function H = clearall(dum, varargin);
% clearall - clear all hoard objects
%    clear(hoard) clear all hoard objects.
%
%    See also hoard, hoard/clear.

clear(dum); % dum serves to make this a hoard method. Thanks & bye.

Shh = get(0,'showhiddenhandles');
set(0,'showhiddenhandles','on');
hmarks = findobj(0,'tag','I_am_a_lonesome_hoard');
set(0,'showhiddenhandles',Shh);
for ii=1:numel(hmarks),
    hfig = get(hmarks(ii), 'parent');
    delete(hfig);
end








