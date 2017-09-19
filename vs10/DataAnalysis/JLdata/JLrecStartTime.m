function T0 = JLrecStartTime(icell_run);
% JLrecStartTime - start time of first recording of given cell
%   T0 = JLrecStartTime(icell_run)
%   in datenum units (i.e. days)
%
%   See also JLrecStartTime

% try cached value
[T0, CFN, CP] = getcache(mfilename, icell_run);
if ~isempty(T0),
    return;
end

% all recs from cell
DB = JLdbase;
DB = DB([DB.icell_run]==icell_run);
DD = JLdatastruct(DB);
T0 = min([DD.abfdatenum]);
putcache(CFN, 100, CP, T0);






