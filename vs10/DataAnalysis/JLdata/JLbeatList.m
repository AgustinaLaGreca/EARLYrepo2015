function S = JLbeatList(ExpID, icell);
% JLBeatList - stripped JLbeat output - all data
%    S = JLBeatList(ExpID);
%    S = JLBeatList;

[ExpID, icell] = arginDefaults('ExpID/icell', []);
% checklist
[L, SL] = JLreadBeats('list');

[S, CFN, Param] = getcache(mfilename, {ExpID icell});
if ~isempty(S),
    return;
end

if ~isempty(ExpID),
    S = JLbeatList;
    ihit = strmatch(upper(ExpID), upper({S.ExpID}));
    S = S(ihit);
    if ~isempty(icell),
        S = S([S.icell]==icell);
    end
else, % all experiments
    % simply load all files previously saved by JLreadBeats
    DD = dir(fullfile(processed_datadir, 'JL', 'JLbeat', 'RG*.mat'));
    S = [];
    for ii=1:numel(DD),
        disp(['---' DD(ii).name '---'])
        qq = load(fullfile(processed_datadir, 'JL', 'JLbeat', DD(ii).name));
        FNS = fieldnames(qq);
        for jj=1:numel(FNS),
            fns = FNS{jj}
            s = local_strip(qq.(fns), fns);
            S = [S; s(:)];
        end
    end
    % check whether all cells are represented
    for ii=1:numel(SL),
        ihit = intersect(strmatch(SL(ii).ExpID, upper({S.ExpID})), find(SL(ii).icell==[S.icell]));
        if isempty(ihit),
            error(['Lost data for ' SL(ii).ExpID ' cell ' num2str(SL(ii).icell) '.']);
        end
    end
end

putcache(CFN, 100, Param ,S);

%=============================
function s = local_strip(s, Filename);
FNS  ={'Stim1'    'Stim2'    'R0'  'Y1'    'Y2'    'Y3'    'Yc' ...
    'Sb'  'Fdpb'  'SPTraw'  'SPT'  'Tsnip'  'Snip'  'APslope'};
for ii=1:numel(FNS),
    [s.(FNS{ii})] = deal([]);
end
[s.seriesID] = deal(Filename);




