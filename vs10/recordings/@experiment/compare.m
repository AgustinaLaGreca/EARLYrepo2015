function T=compare(E,F);
% experiment/compare - compare two experiment definitions
%   compare(E1, E2) lists those paremeters of E2 in that consitute a
%   difference in the definitions of E1 and E2. Differences in the Status
%   field do not count as definition differences.
%
%   L = compare(E1, E2) suppresses listing to the screen and returns a char
%   array T instead.
%
%   See also Experiment/defbackup.

%         ID: [1x1 struct]
%         Audio: [1x1 struct]
%       Electra: [1x1 struct]
%     Recording: [1x1 struct]
%        Status: [1x1 struct]
%       Version: 1

T = '';

FNS = fieldnames(experiment);
for ii=1:numel(FNS),
    fn = FNS{ii};
    if ismember(fn,{'Status', 'Version'}), continue; end % skip these
    [dum, Str] = structcompare(E.(fn), F.(fn));
    if ~isempty(Str),
        T = strvcat(T, ['--------------' fn '------------'], disp(Str));
    end
end
if ~isempty(T),
    T = strvcat(T,'-------------------------------------');
end

if nargout<1,
    disp(T);
    clear T;
end








