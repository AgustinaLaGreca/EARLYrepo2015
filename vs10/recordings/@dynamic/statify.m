function E=statify(E);
% Dynamic/statify - make dynamic object static
%   D = Dataset(D) makes D static, i.e. sets D.isstatic to TRUE.
%
%   See Dynamic.

E.isstatic=1;
upload(E,E.access,'grace');





