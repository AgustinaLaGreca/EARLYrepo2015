function T = JLbeatPool(S);
% JLbeatPool - pool JLbeat analyses
%    S = JLbeatPool([S1 S2 ...]);

Freq1 = unique([S.Freq1]);
Freq2 = unique([S.Freq2]);
SPL1 = unique([S.SPL1]);
SPL2 = unique([S.SPL2]);
if any([numel(Freq1) numel(Freq2) numel(SPL1) numel(SPL2)]>1),
    error('Cannot pool data obtained w different stimuli');
end
FNS = fieldnames(S(1));
T = S(1);
for ii=1:numel(FNS),
    fn = FNS{ii};
    if isequal('RecID', fn);
        T.(fn) = [S.(fn)];
    elseif isequal('ABFname', fn);
        T.(fn) = [S.(fn)];
    elseif  isequal('SPT', fn);
        T.(fn) = cat(1, S.(fn));
    elseif  isequal('SPTraw', fn);
        T.(fn) = cat(1, S.(fn));
    elseif  isequal('Snip', fn);
        T.(fn) = cat(2, S.(fn));
    elseif  isequal('AC', fn);
        T.(fn) = mean([S.(fn)]);
    elseif  isequal('Yoffset', fn);
        T.(fn) = mean([S.(fn)]);
    elseif isnumeric(T.(fn)) && numel(T.(fn)>1000),
        qq = cat(2, S.(fn));
        T.(fn) = mean(qq,2);
    end
end




