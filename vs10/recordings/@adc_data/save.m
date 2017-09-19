function AD = save(AD, DS, FieldName);
% adc_data/save - save ADC data
%   save(AD, DS, FieldName) saves ADC data in chunks and replaced the
%   data in AD by file information on these data. Nothing is done if
%   AD.NsamChunk is infinite, i.e., chunk-saving is disabled.
%   The filename is constructed using dataset/fileprefix.
%
%   See also Dataset/save, dataset/fileprefix, grabbeddata/save.

NumLen = 3; % zero padding constant for file anmes (see below)

Sam = AD.Data.Samples;
if isinf(Sam.NsamChunk), return; end % leave data in place

if ~isempty(Sam.data), % save remaining samples, if any
    % bookkeeping & scaling
    ichunk = Sam.NchunkSaved+1;
    Xmax = max(abs(Sam.data));
    if Xmax==0, Xmax=1; end; % avoid dividing by zero
    Sam.ScaleFactor = [Sam.ScaleFactor, Xmax/(2^31-1)];
    Sam.NchunkSaved = ichunk;
    Filename = [Sam.Filename '_' num2str(ichunk) '.bin'];
    Chunk = int32(Sam.data/Sam.ScaleFactor(end));
    Sam.NsamSaved = [Sam.NsamSaved numel(Chunk)];
    % save samples
    Filename = [Sam.Filename '_' num2str(ichunk) '.bin'];
    fid = fopen(Filename, 'wb');
    fwrite(fid, Chunk, 'int32');
    fclose(fid);
    Sam.data = [];
    AD.Data.Samples = Sam;
end

oldPrefix = [AD.Data.Samples.Filename '_']; % hoard filename - temporary
newPrefix = [fileprefix(DS), '_' FieldName, '_'];
for ichunk=1:AD.Data.Samples.NchunkSaved,
    oldFilename = [oldPrefix, num2str(ichunk) '.bin'];
    newFilename = [newPrefix, zeropadint2str(ichunk,NumLen) '.bin'];
    movefile(oldFilename, newFilename);
end

[Ddir,  AD.Data.Samples.Filename] = fileparts(newPrefix);
AD.Data.Samples.Dir = strrep(lower(Ddir), lower(datadir), ''); % only save subdir relative to datadir
AD.Data.Samples.NumLen = NumLen; % needed to access individual binary files








