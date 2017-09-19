function S = CalcTHR(DS, varargin)
%CALCTHR   calculate threshold curve
%   S = CALCRATE(DS) calculates threshold curve for dataset DS. A structure
%   S is returned containing the calculated data.
%
%   Optional properties and their values can be given as a comma-separated
%   list. To view list of all possible properties and their default value, 
%   use 'factory' as only input argument.

%B. Van de Sande 08-07-2005

%------------------------Default parameters-----------------------
DefParam.isubseqs = 'all'; %The subsequences included in the analysis.
                           %By default all subsequences are included ('all') ...
DefParam.thr      = 10;    %Threshold in dB for the bandwidth and Q-factor ...
DefParam.runav    = 0;     %Number of datapoints used in smoothing the
                           %threshold curve. Smoothing changes the extracted
                           %calculation parameters ...

%-----------------------------------------------------------------
%Checking input parameters ...
if (nargin < 1)
    error('Wrong number of input arguments.'); 
end
if (nargin == 1) & ischar(DS) & strcmpi(DS, 'factory')
    disp('Properties and their default values:');
    disp(DefParam);
    return
end
if ~isa(DS, 'dataset') | ~any(strcmpi(DS.stimtype , {'thr', 'th'}))
    error('First argument should be dataset containing threshold curve.'); 
end

%Checking additional property-list ...
Param = checkproplist(DefParam, varargin{:});
Param = ExpandParam(DS, Param);

%Actual calculation ...
Curve = CalcTHRCurve(DS, Param);    

%Reorganizing output and making sure all fieldnames are lowercase ...
S = struct('calcfunc', mfilename, 'dsinfo', dataset, 'param', Param, ...
    'curve', lowerFields(Curve));

%-------------------------local functions-------------------------
function Param = ExpandParam(DS, Param)

if ~isnumeric(Param.thr) | (length(Param.thr) ~= 1) | (Param.thr <= 0)
    error('Invalid value for property thr.'); 
end
if ~isnumeric(Param.runav) | (length(Param.runav) ~= 1) | (Param.runav < 0) | (Param.runav > length(Param.isubseqs)) 
    error('Invalid value for property runav.'); 
end

%-----------------------------------------------------------------
function Curve = CalcTHRCurve(DS, Param)

%Calculating CF and threshold at CF ...
Freq = DS.Stim.Fcar;
Thr  = DS.Data.Thr;
if elem(size(Freq),2)
   Freq = Freq(:,1); 
end
%Extracting CF and minimum threshold ...
[minThr, idx] = min(Thr);
CF = Freq(idx);

%Extracting Spontaneous Rate ...
if isfield(DS.Data,'SR')
    SR = DS.Data.SR.SR;
else
    SR = -1;
end

%Calculating Q and BW ...
BWThr = minThr + Param.thr;
warning('off','MATLAB:interp1:NaNinY'); %reduce output spam
BWf   = cintersect(Freq, -Thr, -BWThr);
warning('on','MATLAB:interp1:NaNinY'); 
BW    = abs(diff(BWf));
Q     = CF/BW;

Curve = CollectInStruct(Freq, Thr, CF, minThr, SR, BWf, BW, Q);

%-----------------------------------------------------------------