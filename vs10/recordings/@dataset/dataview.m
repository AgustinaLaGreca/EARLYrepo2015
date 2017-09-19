function figh = dataview(DS, P, figh);
% Dataset/dataview - view dataset using a dataviewparam object
%    figh = dataview(DS, P) views dataset DS using the "viewer" dataview(P) 
%    and the parameters in dataviewparam object P. Dataview returns a 
%    graphics handle to the figure.
%
%    dataview(DS, P, figh) uses figure handle figh for the view. By default
%    a new figure is opened.
%
%    See dataviewparam, dataset/listdataviewer, dataset/isdataviewer, 
%    dataset/dotraster, dataset/enableparamedit.

if nargin<3, figh=figure; end % default: open new figure

if ~isa(P, 'dataviewparam'),
    error('Second input argument must be a Dataviewparam object.')
end
% delegate to the dataviewer
feval(dataviewer(P), DS, figh, P);








