function [L, Cstr] = JLdatasearch(ExpID, iCell);
%   JLdatasearch - list all recordings of given cell
%   [L, Cstr] = JLdatasearch('Foo', iCell) returns list of all recordings 
%   from the cell # iCell of experiment Foo. The two output formats are
%       L: char array, each row is a dataset descr
%    Cstr: cell string, same organization.
%
%   See also JLdataList

JLp = JLpool;
imatch = strmatch(lower(ExpID), lower({JLp.SGSRExpID}), 'exact');
JLp = JLp(imatch);
imatch = strmatch([num2str(iCell) '-'], {JLp.SGSRidentifier});
JLp = JLp(imatch);
Nds = numel(JLp);
Cstr = '';
for ii=1:Nds,
    RecID = JLp(ii).SGSRidentifier;
    try,
        DS = sgsrdataset(ExpID, RecID);
        cs = [[ExpID ' '] DSinfoString(DS)];
        cs = cs([1 3:8]);
    catch,
        cs = {[ExpID ' ']  ['<' RecID '>']  '****'   '****'   '****'   '****'   '****'   };
    end
    Cstr = [Cstr; cs];
end
L = '';
for icol = 1:size(Cstr,2),
    L = [L, strvcat(Cstr(:, icol))];
end

% {JLp.SGSRExpID};
% {JLp.SGSRidentifier};
%SGSRdataset(JLp(1).SGSRExpID, JLp(1).SGSRidentifier)









