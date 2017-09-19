function E=finalize(E);
% Dynamic/finalize - get final value of dynamic object
%   D = finalize(D) retrieves the current common value of D and removes it
%   from the cental storage.
%
%   See Dynamic.

E = download(E);
offload(E);
E.isstatic=1;





