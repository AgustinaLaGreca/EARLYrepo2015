function BW = GetNoiseBW(ds, varargin)
%GETNOISEBW return noise bandwidth for dataset
%   BW = GETNOISEBW(ds) returns the bandwidth of the noise used in the
%   supplied dataset as a two-element vector. If the requested playback
%   channel for the dataset wasn't used to administrate a noise token,
%   then NaN is returned.
%
%   Optional properties and their values can be given as a comma-separated
%   list. To view list of all possible properties and their default value, 
%   use 'factory' as only input argument.

%B. Van de Sande 19-07-2004

%-------------------------------default parameters---------------------------
DefParam.channel = 1;   %'master', 'slave', 'left', 'right', 1 or 2 ...
DefParam.numprec = 5e1; %numerical precision for noise bandwidth ...

%----------------------------------main program------------------------------
%Evaluate input arguments ...
if (nargin == 1) & ischar(ds) & strcmpi(ds, 'factory'),
    disp('Properties and their factory defaults:');
    disp(DefParam);
    return;
elseif ~isa(ds, 'dataset'), error('First argument should be dataset.'); end

Param = checkproplist(DefParam, varargin{:});
Param = CheckParam(Param);


BW = GetNoiseBW_EARLY(ds, Param);

%------------------------------local functions-------------------------------
function Param = CheckParam(Param)

if ischar(Param.channel) & any(strncmpi(Param.channel, {'m', 'l'}, 1)), Param.channel = 1;
elseif ischar(Param.channel) & any(strncmpi(Param.channel, {'s', 'r'}, 1)), Param.channel = 2;   
elseif ~isnumeric(Param.channel) | ~any(Param.channel == [1, 2]),
    error('Property channel must be ''master'', ''slave'', ''left'', ''right'', 1 or 2.'); 
end
if ~isnumeric(Param.numprec) | (length(Param.numprec) ~= 1) | (Param.numprec <= 0),
    error('Invalid value for property numprec.');
end

%----------------------------------------------------------------------------
function BW = GetNoiseBW4EDF(ds, Param)

DSSidx = Param.channel;
if (ds.dssnr == 0), error('Dataset doesn''t contain spiketime data.'); end
if (DSSidx > ds.dssnr), error('Slave DSS not used for EDF dataset.'); end

if any(strcmpi(ds.DSS(DSSidx).Mode, {'gws', 'gwr', 'gam'})),
    %Loading the general waveform dataset ...
    FileName = cellstr(ds.GWParam.FileName);
    [dummy, dummy, GWFile] = unravelVMSPath(FileName{DSSidx});
    GWID = cellstr(ds.GWParam.ID);
    GWID = GWID{DSSidx};
    try, GEWAVds = dataset(GWFile, GWID);
    catch, error(sprintf('The general waveform dataset %s <%s> cannot be found.', GWFile, GWID)); end    
    
    %Using EVALGEWAV.M to extract the bandwidth ...
    BW = EvalGEWAV(GEWAVds, 'numprec', Param.numprec, 'verbose', 'no', 'plot', 'no');
else, BW = NaN; end

%----------------------------------------------------------------------------
function BW = GetNoiseBW_EARLY(ds, Param)

if (lower(ds.Stim.DAC(1)) == 'b')
    Nchan = 2;
else 
    Nchan = 1;
end

if isnan(Nchan), error('Dataset doesn''t contain spiketime data.'); end
if (Param.channel > Nchan), error('Second playback channel not used for dataset.'); end

BW = [ds.Stim.LowFreq ds.Stim.HighFreq];

if (size(BW, 1) > 1), BW = BW(Param.channel, :); end

%----------------------------------------------------------------------------