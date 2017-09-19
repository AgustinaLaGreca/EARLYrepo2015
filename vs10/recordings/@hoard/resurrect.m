function H = resurrect(dum,UD);
% hoard/resurrect - reconstruct hoard object from its userdata.
%    resurrect(hoard,UD) retuirns a hoard object with userdata UD.
%
%    See also hoard/getuserdata, hoard/get.

clear(dum); % dum serves to make this a hoard method. Thanks & bye.
H = hoard;
set(H.h, 'userdata', UD);
for ichunk=1:UD.NchunkSaved,
    FN = [UD.Filename '_' num2str(ichunk) '.bin'];
    if ~exist(FN, 'file'),
        error(['Missing hoard file ''' FN '''.']);
    end
end






