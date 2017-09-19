function iv = isdataviewer(dum, M);
% Dataset/isdataviewer - True for dataviewer 
%    isdataviewer(dataset(), @Foo) returns Tru when Foo is a valid 
%    dataviewer. A valid dataviewer is a Dataset method that meets the 
%    requirements decribed in the help text of dataviewparam. 
%    Isdataviewer simply checks whether Foo is in the list of dataviewers 
%    previously registered using Listdataviewer.
%
%   See dataviewparam, dataset/listdataviewer, dataset/dotraster.

L = lower(listdataviewer(dum));
M = lower(char(M)); % use name, not function handle
iv = ismember(M, L);




