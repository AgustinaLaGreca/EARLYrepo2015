function clear(H,flag);
% hoard/clear - clear data of hoard object
%    clear(H) deallocates the memory allocated by H. After clearing H, you 
%    cannot use H anymore. Note that this is 
%    more rigorous than set(H,[]), which just empties H's contents.
%
%    clear(H, 'full') also removes any data on disk that belong to H. This
%    prevents any subsequent resurrection of H.
%
%    See also hoard/set, hoard/cat, hoard/resurrect.

if nargin<2, flag=''; end

UD = getuserdata(H);
delete(H.figh);
EXC = ['builtin(''clear'', ''' inputname(1) ''');'];
evalin('caller', EXC);

if isequal('full', flag),
    delete([UD.Filename '_*.bin']);
end






