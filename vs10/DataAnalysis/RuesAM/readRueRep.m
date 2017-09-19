function [Y, Ncond, dt] = readRueRep(FN, irep, Nd, doDetrend);
% readRueRep - read single rep of Rue data, all conditions concatenated
%   [Y, Ncond, dt] = readRueRep(FN, irep, Mdown, doDetrend);
%    Inputs 
%          FN: filename, e.g. p81213Ac or p81213A (see below).
%        irep: rep number betwen 1 and Nrep
%       Mdown: downsample factor.
%   doDetrend: specify whether piecewise detrending is applied. A "piece"
%              is a region of the recording having the same SPL.
%    Outputs 
%           Y: matrix whose columns are the conditions
%       Ncond: # conditions
%          dt: sample period in ms (taking downsampling into account).
%   If filename FN contains a trailing 'i' of 'c', the contra and ipsi
%   trace with 50-ms baseline removed are used. If not, the complete ipsi & 
%   contra recordings are read and interleaved.
%
%   See also ReadRueCond, RueNrep.

if nargin<3, Nd=1; end
if nargin<4, doDetrend=0; end

if isequal('i', lower(FN(end))) || isequal('c', lower(FN(end))),
    doInterleave = 0;
else,    
    doInterleave = 1;
end
NsamCond = 12500; % # samples per condition (~500 ms)
if isequal('Cel',CompuName),
    error('WHERE?');
else,
    % Ddir = 'D:\Data\RueData\IBW\superRecordings';
    Ddir = 'D:\Data\RueData\IBW2';
end
if ~isdivider(Nd, NsamCond),
    error('Downsampling factor does not divide the number of samples per condition.');
end

Nfreq = 29;
NSPL = 9;
Ncond = Nfreq*NSPL;
Nrep = RueNrep(FN);
if irep>Nrep,     
    error('irep exceeds #reps'); 
end

% read both contra & ipsi recs, restore their original order & detrend
FFNc = fullfile(Ddir, [FN 'c_' num2str(irep-1) '_.ibw']);
FFNi = fullfile(Ddir, [FN 'i_' num2str(irep-1) '_.ibw']);
D = IBWread(FFNc);
Y = reshape(D.y, NsamCond, Ncond);
D = IBWread(FFNi);
meanY = mean(Y(:));
Y = [Y, reshape(D.y, NsamCond, Ncond)];
Y = Y(:, VectorZip(1:Ncond, Ncond+(1:Ncond))); % interleave c-i-c-i...
Y = decimate(Y(:), Nd);
if doDetrend,
    idt = 0:(NsamCond*NSPL*2/Nd):numel(Y);
    idt(1) = 1;
    Y = detrend(Y(:), 'linear', idt);
end
Y = reshape(Y, [NsamCond/Nd 2*Ncond]);
% restore original mean value
Y = Y + meanY - mean(Y(:));
dt = Nd*D.dx;

% if only contra or ipsi is needed, extract it
if ~doInterleave,
    contra = isequal('c', lower(FN(end)));
    NsamCond = 11250; % # samples per condition (~450 ms)
    if contra, % odd columns
        Y = Y(:,1:2:end);
    else, % ipsi: even columns
        Y = Y(:,2:2:end);
    end
end




