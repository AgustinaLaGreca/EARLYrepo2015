function P=StimGuiViewPulldown(src, dum, ID, varargin);
% StimGuiViewPulldown - creator & callback of View Pulldown menu of StimGUI
%    P=StimGuiViewPulldown returns the PulldownMenu object to be
%    included in a StimGUI. It contains these menu items: 
%      Plot waveforms 
%
%    StimGuiViewPulldown also serves as the callback function of the very
%    menu items created by it.
%
%    See also StimGUI, StimPlot, StimGuiEditPulldown, StimGuiFilePulldown.

if nargin<1, % create call
    callme = fhandle(mfilename); % function handle to this function 
    P=pulldownmenu('View','&View');
    P=additem(P,'&Plot waveforms', {callme 'plot'}, 'accelerator', 'L');
else, % callback
    figh = parentfigh(src); % GUI figure handle
    if isequal('plot', ID), 
        % ==========prompt user for param file & save=========
        StimPlot(figh);
    else, 
        error('Invalid callback.');
    end
end
