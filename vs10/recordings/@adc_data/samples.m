function [Y, dt, t0, Dtp] = samples(D, i0, i1, dataSet);
% adc_data/samples - samples of ADC data object.
%   Samples(D) returns the samples adc_data object D in a column array.
%
%   Samples(D, i0, i1) returns a subselection of the samples starting at
%   i0 and ending at i1 (Matlab base-1 convention of indices). i1=inf means
%   up to the last sample. When using this form of Samples, only those 
%   binary files are opened that contain the desired selection. This is 
%   useful to prevent memory problems in the case of very long recordings.
%
%   [Y, dt, t0, DataType] = Samples(D, ...)  also returns the sample 
%   period dt [ms], the timing offset t0[ms] and the DataType char string 
%   (e.g., 'Microphone1'). The timing offset is minus the total delay of 
%   the D/A and A/D conversion. 
%
%   See also adc_data/Fsam, adc_data/Nsam, adc_data/dt, adc_data/Dtype.

if nargin<2, i0=1; end % default: start at first sample
if nargin<3, i1=inf; end % default: end @ last sample

S = D.Data.Samples;
data = S.data;
dt = 1e3/Fsam(D);
t0 = -D.Data.TimingCalib.Lag_ms;
Dtp = dataType(D);

try
    Ddir = datadir(D); % full dir name 
catch
    if i0 > length(data)
        Y=data;
    else
        if ~isinf(i1)
            Y=data(i0:i1);
        else
            Y=data(i0:end);
        end
    end
    return;
end
if isempty(Ddir),
    [dum, Ddir]  = fileparts(S.Dir);
    error(['Cannot find data folder ''' upper(Ddir) '''.']);
end
% determine which fbin files to read from
if isinf(i1),
    i1 = Nsam(D);
end
Y = [];
if i1>Nsam(D),
    error('i1 exceeds # samples saved.');
end
if i0>i1,
    return;
end
if i0<0,
    error('i0 must be positive integer.');
end
% start and end chunk
ich0 = 1+floor((i0-1)/S.NsamChunk); % chunk # of first sample
ipos0 = 1+rem(i0-1,S.NsamChunk); % start index within start chunk 

ich1 = 1+floor((i1-1)/S.NsamChunk);  % chunk # of last sample
ipos1 = 1+rem(i1-1,S.NsamChunk); % stop index within last chunk 

for ich=ich0:ich1,
    FFN = fullfile(Ddir, [S.Filename, zeropadint2str(ich, S.NumLen), '.bin']);
    if ich==ich0, isam0 = ipos0;
    else, isam0=1;
    end
    if ich==ich1, isam1 = ipos1;
    else, isam1=S.NsamChunk;
    end
    [fid, Mess] = fopen(FFN, 'rb');
    error(Mess);
    byte_offset = 4*(isam0-1); % 4 bytes per sample
    Nsam_read = isam1-isam0+1;
    fseek(fid, byte_offset, 'bof');
    y = fread(fid, Nsam_read, 'int32');
    if ~isequal(numel(y), Nsam_read),
        error(['Error while reading binary data file ''' FFN '''']);
    end
    fclose(fid);
    Y = [Y; S.ScaleFactor(ich)*y];
end


