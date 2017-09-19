function c = stimGUIcolor(StimName);
% stimGUIcolor - color of stimulus GUI
%    stimGUIcolor('Foo') returns background color of stimulus GUI for
%    stimulus type Foo.

switch upper(StimName),
    case 'FS', c = [1 0 0];
    case 'NPHI', c = [0 0 1];
    otherwise, c = GUIsettings('GUIfigure')>'Color';
end

