function typelog(E, Ename);
% experiment/typelog - display text from experiment log file
%    typelog(E) displays the text of E's log file. 
%
%    typelog(experiment(), 'foo') displays the text of Exeriment named Foo.
%
%    See also experiment/addtolog, experiment/insertnote,
%    experiment/status, experiment/stimlist, textWrite.

if isvoid(E)
    if nargin<2, 
        error('Cannot type log of void experiment.');
    else,
        E = find(E,Ename);
    end
end
type([filename(E) '.log']);












