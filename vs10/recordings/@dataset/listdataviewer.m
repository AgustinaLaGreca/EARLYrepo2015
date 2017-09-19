function L = listdataviewer(dum, kw, M)
% Dataset/listdataviewer - list of dataviewers for dataset analysis
%    ListDataviewer(dataset()) returns the names of dataset methods that
%    are currently registered as "dataviewers". A dataviewer is a dataset
%    method that meets the syntax requirements described in the help text
%    of dataviewparam.
%
%    ListDataviewer(dataset(), 'clear')  empties the registered list of
%    dataviewers.
%
%    ListDataviewer(dataset(), 'add', @Foo) adds Foo to the list. To test
%    whether Foo is a worthy dataviewer candidate, a call to dataviewparam
%    is employed. 
%
%    See also dataviewparam, dataset/isdataviewer, dataset/dotraster.

if nargin<2, kw = 'get'; end % query

% where on disk is the list kept? 
CacheDir = fileparts(which([class(dum), '/' mfilename]));
CFN = fullfile(CacheDir, 'DataviewerList.EarlySetup');
Cparam = 1; % fake cache parameter.

switch kw,
    case 'get', % query
        L = getcache(CFN,Cparam);
    case 'clear';
        L = {};
        putcache(CFN, 1, Cparam, L);
    case 'add',
        try,
            dataviewparam(M);
        catch,
            error(['Unable to create dataviewparam object from dataset method ''' char(M) '''.']);    
        end
        L = getcache(CFN,Cparam);
        if isempty(L),
            L = {char(M)};
        else,
            L = [L char(M)];
            L = unique(L);
        end
        putcache(CFN, 1, Cparam, L);
end




















