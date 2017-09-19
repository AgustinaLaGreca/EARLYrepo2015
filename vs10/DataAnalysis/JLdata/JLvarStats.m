function S = JLvarStats(ExpID, iCell);
% JLvarStats - stats from JLbeatVar
%   S = JLvarStats('RG09178', 1)  % single cell
%   S = JLvarStats;  % all cells

if nargin<1, % all
    [dum L] = JLreadBeats('list');
    S = [];
    for ii=1:numel(L),
        disp(L(ii))
        S = [S JLvarStats(L(ii).ExpID ,L(ii).icell)];
    end
    return;
end

% ---single cell from here---
% cache
[S, CFN, CP] = getcache(mfilename, {ExpID, iCell});
if ~isempty(S), 
    return; 
end

S = [];
JLreadBeats(ExpID, iCell);
DS = who('Jb*');
for ids = 1:numel(DS),
    Jb = eval(DS{ids});
    for icond=1:numel(Jb),
        s = JLbeatVar(Jb(icond),0);
        s = local_strip(s,1);
        S = [S s];
    end
end
putcache(CFN, 200, CP, S);

%=========================
function S = local_strip(S, Rig);
Rig = arginDefaults('Rig',0);
for ii=1:numel(S),
    s = S(ii);
    s.bin_MeanWave = [];
    s.bin_VarWave = [];
    if Rig,
        s.ipsi_MeanWave = [];
        s.ipsi_VarWave = [];
        s.contra_MeanWave = [];
        s.contra_VarWave = [];
        s.sil_MeanWave = [];
        s.sil_VarWave = [];
    end
    S(ii) = s;
end


















