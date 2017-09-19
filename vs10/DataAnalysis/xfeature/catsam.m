function [sam, T] = catsam(D, isweep, ichan);
% catsam - concatenate all samples of TK ABF file
%    [sam, Time] = catsam(D, isweep, ichan)
%    Default all sweeps, ichan=1.
%    dt in ms.
%
%    [sam, Time] = catsam({'RG09178' '1-19-'}, isweep, ichan)

if nargin<2, isweep=1:D.Nsweep; end
if nargin<3, ichan=1; end

if iscell(D),
    D = readTKABF(D{:}),
elseif isstruct(D) && ~isfield(D,'AD'),
    D = readTKABF(D);
end

sam = cat(1,D.AD(isweep,ichan).samples); 
T = timeaxis(sam, D.dt_ms);
