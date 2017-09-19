function [D, DS, E] = clampfit(ExpID, RecID, icond, varargin);
% clampfit - open ABF data in Clampfit
%   clampfit(ExpID, RecID, icond) opens the ABF file behind the dataset
%   specified by (ExpID, RecID, icond). See
%
%   clampfit(S) opens the ABF file behind the dataset specified by 
%   S.ExpID, S.RecID, S.icond. S can be the ABF dataset itself (see
%   readTKabf) or the output of an analysis function like JLbeat.
%
%   clampfit() uses the most recently loaded data.
%
%   [D, DS, E] = clampfit(...) returns ABF dataset D, SGSR dataset DS and
%   pool element E as returned by readTKabf.
%
%   See also readTKabf, lastABF.

if nargin<1,
    %qq=myflag('LastRecordingAnalyzed');
    qq=MyFlag('LastRecordingLoaded');
    [ExpID, RecID, icond] = deal(qq{:});
end
if nargin==1 && isstruct(ExpID), % a struct, hopefully containing rec info
    [ExpID, RecID, icond] = deal(ExpID.ExpID, ExpID.RecID, ExpID.icond);
elseif nargin==1 && iscell(ExpID), % distribute
    [ExpID, RecID, icond] = deal(ExpID{:});
elseif nargin==1 && isnumeric(ExpID), % unique cell number
    ExpID = JLdatastruct(ExpID);
    [ExpID, RecID, icond] = deal(ExpID.ExpID, ExpID.RecID, ExpID.icond);
end
[D, DS, L]=readTKABF(ExpID, RecID, icond);

FN = fullfile(L.Dir, L.Recs{icond});
if ~isequal('Kulak', CompuName),
    %D = readABFdata(FN);
    fh = figure;
    set(fh,'units', 'normalized', 'position', [0.261 0.13 0.706 0.777])
    switch D.Nchan,
        case 4, ABFplot(D,0,[1 3 4]); % omit void chan #2
        otherwise 2, ABFplot(D); % plot all channels
    end
else,
    winopen(FN);
end







