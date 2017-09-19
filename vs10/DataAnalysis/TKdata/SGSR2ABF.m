function [ABFdir, ABFfiles, FullABFdir, FullABFname] = ABF_fileInfo(S, ExpID, RecID, irep);
% ABF_fileInfo - return ABF datadir and filenames corresponding to single SGSR recording
%   [ABFdir, ABFfiles, FullABFdir] = SGSR2ABF(Sall, ExpID, RecID);
%   returns directory and all ABF filenames, sorted according to date.
%
%   [ABFdir, ABFfile, FullABFdir, FullABFname] = SGSR2ABF(Sall, ExpID, RecID, irep);
%   only returns the specified ABF filename. irep=0 is the first, 
%   "spontaneous" recording.
%
%   Examples
%      Sall = load('D:\USR\Marcel\sortTKdata\all.SGSR_ABFs', '-mat'); 
%      Sall=Sall.S;
%      [ABFdir, ABFfiles, FullABFdir] = SGSR2ABF(Sall, 'RG09110', '3-1-FS')
%   ABFdir =
%   050509\pen3u7\recording39
%   ABFfiles = 
%       '050509Pen1_0860.abf'    '050509Pen1_0861.abf'    ...
%   FullABFdir =
%       D:\USR\Thomas\data\050509\pen3u7\recording39
%  
%     [ABFdir, ABFfiles, FullABFdir, FN] = SGSR2ABF(Sall, 'RG09110', '3-1-FS', 0)
%   ...
%   FN =
%       D:\USR\Thomas\data\050509\pen3u7\recording39\050509Pen1_0860.abf

% select exp
iexp = strmatch(ExpID, {S.SGSRExpID}, 'exact');
S = S(iexp);
% select rec
irec = strmatch(upper(RecID), {S.SGSRidentifier});
if isempty(irec),
    error(['No ABF recordings matching ''' ExpID, '/' RecID '''.']);
elseif numel(irec)>1,
    error(['Multiple ABF recordings matching ''' ExpID, '/' RecID '''.']);
end
S = S(irec);
ABFdir = fullfile(S.ExpID, S.Pen, S.recID);
FullABFdir = S.Dir;
qq = dir(fullfile(FullABFdir, '*.abf'));
[dum, isort] = sort([qq.datenum]);
ABFfiles = {qq(isort).name};

if nargin>3,
    ABFfiles = ABFfiles{1+irep} % see help text
    FullABFname = fullfile(FullABFdir, ABFfiles);
end







