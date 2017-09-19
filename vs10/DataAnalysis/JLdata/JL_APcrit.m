function APcrit = JL_APcrit(ExpID, iCell);
%  JL_APcrit - find criterium for action potentials of given JL cell.
%     JL_APcrit(ExpID, iCell)

if isstruct(ExpID), % output struct from JLbeat or JLstats
    [ExpID, iCell] = deal(ExpID.ExpID, ExpID.icell);
end

qq = dir(which('JLreadBeats'));
CP = {qq.date ExpID, iCell};
CFN = mfilename;
APcrit = getcache(CFN, CP);
if ~isempty(APcrit),
    return;
end

APcrit = nan;
tag = upper(['local_' ExpID '_' num2str(iCell)]);
LL = mytextread(which('JLreadBeats'));
for ii=1:numel(LL),
    if ~isempty(strfind(upper(LL{ii}), tag)),
        for jj=ii+(1:1000),
            if ~isempty(strfind(upper(LL{ii}), 'function')),
                error('APcrit not found in JLreadBeats.m');
            elseif ~isempty(strfind(LL{jj}(1:min(end,10)), 'APcrit')); 
                % should be a line like: APcrit = [-10 -1 1.4 -7]; Nav = 1; % mediocre thr definition
                eval(LL{jj}); 
                putcache(CFN, 1, CP, APcrit);
                return;
            end
        end
    end
end








