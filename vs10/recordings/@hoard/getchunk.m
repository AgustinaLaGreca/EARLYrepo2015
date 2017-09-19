function C=get(H, ichunk);
% hoard/get - get data of hoard object
%    get(H) returns the data of hoard object H 
%
%    See also hoard/set, hoard/cat, hoard/getchunk.

UD=get(H.h, 'userdata');
if ichunk>UD.NchunkSaved,
    error(['Chunk index exceeds # chunks saved.']);
end
Filename = [UD.Filename '_' num2str(ichunk) '.bin'];
fid = fopen(Filename, 'rb');
C = UD.ScaleFactor(ichunk)*fread(fid, inf, 'int32=>double').';
fclose(fid);




