function C = StimCollection;
% StimCollection - list of all available stimuli
%     StimCollection() returns a cell string C whose elements
%     represent all the stimulus names for which a definition is present in
%     the stimdefPath. A stimulus FOO is represented by a file stimdefFOO
%     in the stimdefPath. Use stimButtonTiling to select the subset that is
%     represented by buttons on the dashboard.
%
%     See also Dashboard, stimButtonTiling, StimdefPath.

C = {};
DD = Words2cell(stimdefPath, pathsep);
for ii=1:numel(DD),
    qq = dir(fullfile(DD{ii}, 'stimdef*.m'));
    names = lower({qq.name});
    names = strrep(names, 'stimdef', '');
    names = strrep(names, '.m', '');
    C = [C upper(names)];
end








