function versiondir = findversiondir(EarlyRoot);
% findversiondir - find dir that contains newest EARLY version 
%   helper function for startup.m

qq = dir([EarlyRoot filesep 'vs*']);
if ~isempty(qq), 
    qq = qq([qq.isdir]); 
end
if isempty(qq),
    error(['Cannot find a candidate EARLY versiondir ''' [EarlyRoot filesep 'vs*'] '''.']);
end

Nam = strvcat({qq.name});
Vs = Nam(:,3:end); % version numbers in char array
[dum imax] = max(str2num(Vs)); % highest version number
versiondir = fullfile(EarlyRoot, qq(imax).name);

