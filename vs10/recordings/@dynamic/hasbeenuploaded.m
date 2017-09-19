function hb=hasbeenuploaded(D);
% Dataset/hasbeenuploaded - true for dynamic uploaded objects
%   hasbeenuploaded(D) tests wether D has been previously uploaded.
%
%   See Dynamic/upload.

hb = ~isempty(D.access);


