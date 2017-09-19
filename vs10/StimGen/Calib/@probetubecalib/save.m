function [okay, FFN] = save(T, FN, Prompt);
% probetubecalib/save - save Probetubecalib object to file.
%    save(T,'Foo') saves T to file Foo. Explicit directory or extension 
%    will overrule the default dir Datadir(T) and default extension 
%    fileExtension(T). If no filename is specified, the user is prompted as 
%    described below.
%
%    [okay, FFN]=save(T, '?', P) prompts the user for a filename using UIPUTFILE with
%    prompt P. The default prompt is 'Specify file name for Transfer data.'
%    Save returns True if T got saved, False otherwise. The latter happens
%    when the user cancels the UIPUTFILE menu. If the user cancels, okay is
%    False and FFN = '', otherwise okay=True and FFN is the full file name.
%
%    save(T, '?Foo', P) prompts the user for a filename applying a
%    filefilter Foo. 
%
%    Examples
%      save(T, 'Foo')
%      save(T, 'Foo.extension')
%      save(T, 'mydir/Foo')
%      save(T, '?*.extension')
%      save(T, '?dir/*.extension')
%
%    See also Probetubecalib/datadir, Probetubecalib/load, UIPUTFILE.

okay = 1;
if nargin<2,
    FN = '?';
end
if nargin<3,
    Prompt = 'Specify file name for Transfer data.';
end
if isequal('?', FN(1)), % prompt for file
    FileFilter = FN(2:end);
    if isempty(FileFilter),
        FileFilter = '*';
    end
    FileFilter = FullFilename(FileFilter, datadir(T), fileExtension(T));
    [FN, PATH]  = uiputfile(FileFilter, Prompt);
    if isequal(0, FN), okay=0; FFN = ''; return; end % user cancelled
    FN = fullfile(PATH, FN);
end
FFN = fullfilename(FN, datadir(T), fileExtension(T));
save(FFN,'T','-mat');





