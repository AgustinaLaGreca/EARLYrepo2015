function [okay, FFN] = save(T, FN, Prompt);
% transfer/save - Transfer object to file.
%    save(T,'Foo') saves T to file Foo.transfer in the calibration
%    directory. Explicit directory or extension will overrule the default
%    dir (calibrationdir) and extension ('transfer').
%
%    save(T, '?', P) prompts the user for a filename using UIPUTFILE with
%    prompt P. The default prompt is 'Specify file name for Transfer data.'
%    save returns True if T got saved, False otherwise. The latter happens
%    when the user cancels the UIPUTFILE menu.
%
%    save(T, '?Foo', P) prompts the user for a filename applying a
%    filefilter Foo. 
%
%    Examples
%      save(T, 'Foo')
%      save(T, 'Foo.extension')
%      save(T, 'mydir/Foon')
%      save(T, '?*.extension')
%      save(T, '?dir/*.extension')
%
%    See also Calibrationdir, Transfer/load, UIPUTFILE.

okay = 1;
if nargin<3,
    Prompt = 'Specify file name for Transfer data.';
end
if isequal('?', FN(1)), % prompt for file
    FileFilter = FN(2:end);
    if isempty(FileFilter),
        FileFilter = '*';
    end
    FileFilter = FullFilename(FileFilter, calibrationdir, '.transfer');
    [FN, PATH]  = uiputfile(FileFilter, Prompt);
    if isequal(0, FN), okay=0; FFN = ''; return; end % user cancelled
    FN = fullfile(PATH, FN);
end
FFN = fullfilename(FN, calibrationdir, '.transfer');
save(FFN,'T','-mat');





