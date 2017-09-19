function [S, ST] = JLgetBeatStruct(S, Stripped);
% JLgetBeatStruct - retrieve full JLbeat output from stripped JLbeatStats output
%    [S, ST] = JLgetBeatStruct(S, Stripped)
%    inputs
%       S: either JLbeat or JLbeatStats output
%   Stripped: logical. If 1, use JLbeatList to get quick srtrpiied version
%          of JLbeat output.
%    output:
%        S: JLbeat output
%       ST: JLbeatStats output
%
%    See also JLfind.

Stripped = arginDefaults('Stripped', 0);
if Stripped,
    qq = JLbeatList();
end


if isfield(S,'Vbin_STD'), % JLbeatStats output; store it
    ST = S;
else,
    ST = [];
end

if ~isfield(S,'R0'), % S is return arg of JLbeatStats. Retrieve original JLbeat output as saved by JLreadBeats.
    SFN = local_filename(S.ExpID, S.icell);
    qq = load([SFN '.mat'],'-mat');
    qq = qq.(['Jb' S.seriesID]);
    icond = [S.icond];
    S = qq(icond);
end

if isempty(ST) && nargout>1,
    ST = JLbeatStats(S);
end
%=====================================================
function [SFN, ID] = local_filename(ExpID, iCell);
% filename for processed data
ID = [upper(ExpID) '_' num2str(iCell)];
SFN = fullfile(processed_datadir, 'JL', 'JLbeat', ID); % filename for processed data





