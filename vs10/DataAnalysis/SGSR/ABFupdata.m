function ABFupdata(Experimenter);
% ABFupdata - update ABF data from Kulak
%   ABFupdata('TK') copies any newly collected data from \\Kulak\data_tk
%   local rawdata\sgsr directory.
%
%   See also SGSRdataset.

if nargin<1, Experimenter='TK'; end

if isequal('KULAK', upper(CompuName)),
    error('Kulak is source of ABF files.');
end

Src = ['\\Kulak\data_' Experimenter '$'];
Dest = TKdatadir(Experimenter);
more off;
qq = dir(Src); 
qq = qq([qq.isdir]);
qq = qq(3:end); % remove .\ and ..\
expDirnames = upper({qq.name})

qq = dir(Dest); 
qq = qq([qq.isdir]);
qq = qq(3:end); % remove .\ and ..\
Copied = upper({qq.name})

tobeCopied = expDirnames(~ismember(expDirnames, Copied));

for ii=1:numel(tobeCopied),
    xd = tobeCopied{ii};
    if ~exist(fullfile(Src, xd ,'logs'), 'dir'),
        disp(['Skipping ''' xd ''': no LOGS sub directory found.']);
    elseif ~isempty(dir(fullfile(Src, xd ,'*.abf'))),
        disp(['Skipping ''' xd ''': not all ABF files have been moved to sub folders.']);
    else,
        [Success,Mess, MessID] = mkdir(Dest, xd);
        error(Mess);
        disp(['Copying ''' xd '''.']);
        [Success,Mess, MessID] = copyfile(fullfile(Src, xd), fullfile(Dest, xd));
        error(Mess);
    end
end







