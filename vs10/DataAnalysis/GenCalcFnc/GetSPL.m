function SPL = GetSPL(ds)
%GETSPL get sound pressure level for dataset
%   SPL = GETSPL(ds) returns the sound pressure level in dB for the dataset ds.
%   The SPL is returned as a matrix with the number of rows defined by the total 
%   number of subsequences and with the same number of columns as active channels 
%   used in the collection of data.

%B. Van de Sande 24-05-2004

if (nargin < 1), error('Wrong number of input parameters.'); end
if ~isa(ds, 'dataset'), error('First argument should be dataset.'); end

SPL = Try2ExtractSPL(ds); 
end %Unknown dataset format ...

%---------------------------------------local functions---------------------------------------
function SPL = Try2ExtractSPL(ds)
%Try to retrieve SPL using the virtual field spl of a dataset

if strcmp(ds.Stim.DAC,'Both')
    NChan   = 2;
else
    NChan   = 1;
end
NSub    = ds.Stim.Presentation.Ncond;


SPL = repmat(NaN, NSub, NChan);

try 
    [NRow, NCol] = size(ds.Stim.SPL);
    SPL(1:NSub, 1:NChan) = repmat(ds.Stim.SPL, (NSub-NRow)+1, (NChan-NCol)+1); 
end 
end

