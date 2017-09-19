function [idx, datenum]=calib_index(E);
% experiment/calib_index - index of most reecent calibration data 
%   [idx, datenum]=calib_index(E) returns the index of most recently saved calibration
%   data for experiment E. Zero is returned if no calibraion data exist.
%
%   See also earcalib/load.

idx = 0; datenum = nan; % defaults
dd = dir(fullfile(folder(E), '*.EarCalib'));
if isempty(dd),
    return;
end
dd = sortAccord(dd,[dd.datenum]);
dd = dd(end); % most recent
[dum, idx] = strtok(dd.name,'_');
idx = strtok(idx(2:end),'.');
idx = str2num(idx);
datenum = dd.datenum;










