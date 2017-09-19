function SGSRupdata;
% SGSRupdata - update exported SGSR data from 3rd floor PC
%   SGSRupdata copies any newly exported SGSR data from 3floor measurement
%   PC (X:\Rdam) to the local rawdata\sgsr directory.
%
%   See also SGSRdataset

Src = 'X:\Rdam';
Dest = SGSRdatadir;

more off;
qq = dir(fullfile(Src, 'exported_*')); qq = qq([qq.isdir]);
expDirnames = upper({qq.name});

qq = dir(fullfile(Dest, 'exported_*')); qq = qq([qq.isdir]);
Copied = upper({qq.name});

tobeCopied = expDirnames(~ismember(expDirnames, Copied));

for ii=1:numel(tobeCopied),
    xd = tobeCopied{ii}
    mkdir(Dest, xd);
    s = fullfile(Src, xd, '*.*');
    d = fullfile(Dest, xd);
    [dum, Mess] = copyfile(s , d);
    warning(Mess);
end




