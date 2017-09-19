function [df, MG, PH, Nsam, Pnoise] = ZWOAEspec(D, idataset);
% ZWOAEspec - spectrum of ZWOAE data
%   [df, Magn, Phase, N, Pnoise] = ZWOAEspec(Gerbil, idataset) returns the
%   spectrum of ZWOAE data (Gerbil, idataset) in terms of the output args
%        df: frequency spacing in kHz
%      Magn: magnitude spectrum in dB SPL
%     Phase: phase in cycles
%         N: # samples in spectrum
%    Pnoise: polynomial coeff to estimate noise floor [dB SPL] @ freq [kHz]
%
%    [...] = ZWOAEspec(D), with D being struct returned by either getZWOAEdata
%    instead of reading the data from file.
%
%    For coustical data that contain a nontrivial 'probename' field, the
%    probe transfer function indicated by that field is used to convert the
%    acoustic spectrum as measured by the probe microphone to the "true"
%    spectrum at the eardrum. This compensation of probe transfer includes
%    both magnitude and phase.
%
%    See also getZWOAEdata, ZWOAEfit.

if ~isstruct(D), % read data from file
    D = getZWOAEdata(D,idataset);
end

if isfield(D,'time') ~= 1,
    %struct D comes from ZWpar
    micgain = D.micgain;
    dBV = D.dBV;
    dBSPL = D.dBSPL;
else
    %struct D comes from getZWOAEdata
    micgain = D.micgain;
    dBV = D.micsensitivity.dBV;
    dBSPL = D.micsensitivity.dBSPL;
end

dt = 1/D.fs; % sample period in ms
Nsam = numel(D.signal);

% spectrum
df = 1/(Nsam*dt); % spectral spacing in kHz
Cspec = fft(D.signal).'; % complex spectrum in row vector
PH = angle(Cspec)/2/pi; % phase in cycles
MG = A2dB(abs(Cspec)*2/Nsam); % magnitude in dB
MG = MG-micgain-dBV+dBSPL;

% noise floor estimation
if (D.periodicity==Nsam) && nargout>4,
    warning('Noise floor cannot be estimated; re-compactify the raw data.');
    Pnoise = nan;
end

if nargout>4,
    Pnoise = fitNoise(MG, Nsam/D.periodicity, df);
end

% take any probe transfer into account for MG and PH, not for NoiseFloor
if strcmpi(D.RecType,'acoustic') && isfield(D,'probename') && ~strcmpi(D.probename,'none'),
    freq = 1e3*(0:Nsam-1).'*df; % freq of spectral components in Hz
    [ProbeGain, ProbePhase] = ProbeLoss(D.probename, freq);
    MG = MG - ProbeGain; % minus sign: we "undo" the probe transfer
    PH = PH - ProbePhase; % ditto
end





