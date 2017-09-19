function [X,Nch]=get(H, flag);
% hoard/get - get data of hoard object
%    get(H) returns the dynamic data of hoard object H. This does NOT
%    include any samples saved to disk.
%
%    get(H, 'full') returns the comple array of samples from H, including
%    those stored on disk. Note that the size of the data may be large,
%    which may cause memory errors.
%    
%
%    See also hoard/set, hoard/cat, hoard/getchunk, hoard/getuserdata.

if nargin<2, flag = ''; end

U=get(H.h, 'userdata');
Nch = U.NsamChunk;
X = U.data;

if isequal('full', lower(flag)),
    for ii=U.NchunkSaved:-1:1,
        X = [getchunk(H,ii) ,X];
    end
elseif isempty(flag), % ok
else, error(['Invalid flag ''' flag '''.']);
end






