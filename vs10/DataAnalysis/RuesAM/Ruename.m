function Ruename;

SourceDir = 'D:\Data\RueData\IBW\Audiograms';
DestDir = 'D:\Data\RueData\IBW\Audiograms\NewNames';
qq = mytextread('D:\Data\RueData\IBW\Audiograms\ruename.txt');
for ii=1:numel(qq),
    W = Words2cell(qq{ii});
    [ofn, nfn] = deal(W{:})
    if ~isequal('-', nfn),
        FileFilter = fullfile(SourceDir, [ofn '*.ibw']);
        fL = dir(FileFilter);
        fL = {fL.name};
        for jj=1:numel(fL),
            oldname = fL{jj};
            newname = strrep(oldname, ofn, nfn);
            oo = fullfile(SourceDir, oldname)
            nn = fullfile(DestDir, newname)
            copyfile(oo,nn);
        end
    end
end











