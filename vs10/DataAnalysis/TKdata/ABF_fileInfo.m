function [ABFdir, ABFfiles, FullABFdir, FullABFname] = ABF_fileInfo(S, icond);
% ABF_fileInfo - return ABF datadir and filenames of single SGSR_ABF database entry
%   [ABFdir, ABFfiles, FullABFdir] = ABF_fileInfo(S);
%   returns directory and all ABF filenames, sorted according to date.
%
%   [ABFdir, ABFfile, FullABFdir, FullABFname] = ABF_fileInfo(S, icond);
%   only returns the specified stimulus condition #icond, as specified
%   in the SGSR stimulus.
%
%   EXAMPLE 1
%      Sall = load('D:\USR\Marcel\sortTKdata\all.SGSR_ABFs', '-mat'); Sall=Sall.S;
%      S = SGSR_ABF_filter(Sall, 'RG09110', '3-1-FS')
%      [ABFdir, ABFfiles, FullABFdir] = ABF_fileInfo(S)
%   ABFdir =
%   050509\pen3u7\recording39
%   ABFfiles = 
%       '050509Pen1_0860.abf'    '050509Pen1_0861.abf'    ...
%   FullABFdir =
%       D:\USR\Thomas\data\050509\pen3u7\recording39
%  
%   EXAMPLE 2
%     [ABFdir, ABFfiles, FullABFdir, FN] = ABF_fileInfo(S, 1)
%   ...
%   FN =
%       D:\USR\Thomas\data\050509\pen3u7\recording39\050509Pen1_0860.abf
%
%   See also readTKABF.

error('obsolete function!!')
ABFdir = 123; % obsolete..
qq = dir(fullfile(TKdatadir, '*.abf'));
[dum, isort] = sort([qq.datenum]);
ABFfiles = {qq(isort).name};
if nargin>1, % only a specific condition
    ABFfiles = ABFfiles{icond}; % see help text
    FullABFname = fullfile(TKdatadir, ABFfiles);
end







