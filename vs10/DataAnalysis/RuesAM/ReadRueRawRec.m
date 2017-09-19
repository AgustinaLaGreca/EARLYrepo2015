function D = ReadRueRawRec(FN, ichan, Vthr, Twin)
% ReadRueRawRec - read raw Rue data and truncate spikes
%    ReadRueRawRec(FN, Vthr, Twin)
%    Vmin is threshold; Twin indicates time window around peak, e.g.,
%    Twin = [-1 1] means from 1 ms before peak to 1 ms after peak.
%
%   See APtruncate.


if nargin<3, Vthr = inf; end % default: no truncation
if nargin<4, Twin = [-1 1]; end
Ddir = 'D:\TEMP\RueCorr\';
FFN = fullfile(Ddir, [FN, '.abf'])
D = readABFdata(FFN);
dt = D.dt_ms;
D.AD = cat(1, D.AD(:,ichan));
D.AD = cat(1, D.AD.samples);

minDistsam = round(0.5*diff(Twin)/dt);
[dum, ipeak] = findpeaks(D.AD+1e-4*std(D.AD)*rand(size(D.AD)), ... % add a little noise to unequalize numbers that are equal by rounding 
    'minpeakheight', Vthr, 'minpeakdistance', minDistsam);
tpeak = dt*(ipeak-1);
D.AD = cutFromWaveform(dt, D.AD, tpeak, Twin);