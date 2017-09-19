function S = stimButtonTiling(C);
% stimButtonTiling - arrangement of stimulus buttons on dashboard
%     C = stimButtonTiling() returns a cell string C whose elements
%     represent the respective rows of stimulus buttons on the dashboard.
%     The names of the buttons are separated by slashes. 
%
%     stimButtonTiling(C) sets the tiling to cell string C. The cell string
%     C will be parsed and the existence of the stimuli is checked before
%     it is accepted. Once accepted, the new tiling persists across Matlab
%     sessions. Note that the "existence" of a stimulus named Foo is
%     equivalent to the presence of a file named stimdefFoo somewhere in
%     the stimdefPath. Use stimCollection to see all avalibale stimuli.
%
%     stimButtonTiling('factory') resets the tiling to the hardcoded
%     factory settings.
%
%     Example
%        stimButtonTiling({'FS/THR/RC'  'NRHO/ARMIN'  'ZW/SUP'})
%
%     See also Dashboard, StimCollection, StimdefPath.


SFN = 'LocalStimuli'; % name of setup file
S_fac = {'FS/RF'  'NPHI/MASK'  'SUP/ZW/BINZW'};
if nargin<1, % get
    S = FromSetupFile(SFN, 'StimButtonTiling', '-default', S_fac);
else, % (re)set
    if isequal('factory',C),
        ToSetupFile(SFN, '-propval', 'StimButtonTiling', S_fac);
        S = S_fac;
        return;
    end
    % real set: first check C
    if ~iscellstr(C),
        error('Input C must be cell array of character strings.');
    end
    Nrow = numel(C);
    for irow=1:Nrow,
        Row = Words2cell(C{irow},'/');
        for ii=1:numel(Row),
            if ~isvarname(Row{ii}),
                error(['Invalid stimulus name ''' Row{ii} '''; not a Matlab identifier (see ISVARNAME).']);
            end
            FN = ['stimdef' Row{ii} '.m']; % filename of stimulus definer 
            if isempty(searchInPath(FN, stimdefPath)),
                error(['Stimulus definer ''' FN ''' not found in stimdefPath.']);
            end
        end
    end
    % C is okay; save it & return it
    ToSetupFile(SFN, '-propval', 'StimButtonTiling', C);
    S = C;
end




