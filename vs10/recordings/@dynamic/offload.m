function E=offload(E);
% Dynamic/offload - offload previosuly uploaded version of Dynamic object
%   offload(D) removes the value of D that was previously stored
%   using upload(D) and returns it (with its address empty).
%   
%   See Dynamic/upload, Dynamic/download.

if isempty(E.access),
    error('Cannot offload object that has never been uploaded.');
end
if nargout>1,
    E = download(E);
end
rm(E.access); 
E.access = [];


