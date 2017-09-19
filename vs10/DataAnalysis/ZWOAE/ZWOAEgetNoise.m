function [MG, freqs] = ZWOAEgetNoise(expID, recID, freqs)
% ZWOAEgetNoise - get noise floor in ZWOAE recordings at given freqs
%
% [MG, freqs] = ZWOAEGETNOISE(expID, recID) returns the estimated noise floor
% (and corresponding freqs) of the specified recording. The returned freqs
% are "wrapped" in the range 0...Nyquist
% 
% To compensate for the probe, it takes any applicable transfer function
% for the recording probe into account. This is when:
%   1) the recording is ACOUSTIC
%   2) the probename is non-trivial (e.g. 'copperboy')
%
% [MG, freqs] = CORRECT_FOR_PROBE(expID, recID, freqs) allows specification of
% frequencies to which to apply the transfer function.
%
% Alternative to numeric IDs for expID and recID, expID may also be a
% struct as obtained via ZWOAEimport.m (or equivalent)
%
% To calculate the noise floor, first the polynomial coeffients Pnf
% as returned by ZWOAEimport.m are evaluated. This gives the noise floor
% BEFORE probe transfer. Second, the probe transfer is applied to the
% magnitudes @ these frequencies.
%
% INPUTS:
%       expID:  numeric ID of the experiment <<OR>> datastruct D
%       recID:  single numeric ID of recording
%       freqs:  (array of) frequencies for which to return the Noise [kHz]
%
% OUTPUTS:
%       MG:     (array of) amplitudes for estimated noise [dB SPL]
%       freqs:  (array of) corresponding freqs [kHz]; in the range 0...Nyquist
%
% Syntax examples:
%       ZWOAEgetNoise(80, 1230);
%       ZWOAEgetNoise(80, 1230,[3 4 5 6]);
%
%       D = ZWOAEimport(80, 1230);
%       ZWOAEgetNoise(D);
%       ZWOAEgetNoise(D,[3 4 5 6]);
%
% See also ZWOAEspec, probeloss, ZWOAEimport
%

if nargin < 3, freqs = []; end

if isstruct(expID), %data given via struct
    if nargin == 2,
        freqs = recID;
    end
    if isfield(expID,'MG'),
        D = expID;
    else,
        error('Specified struct D does not hold recorded data, just params.');
    end
else
    D = ZWOAEimport(expID, recID);
end

if numel(D) > 1, error('Only one recording at a time is allowed.'); end

if isempty(freqs),
    %if omitted, generate frequency axis covering all freqs in spectrum
    halfSam = D.Nsam/2;
    freqs = (0:halfSam-1).'*D.df; % freq of spectral components in kHz
end
%get noise floor from data WITHOUT probe correction
[MG, freqs] = local_NoisefromData(D, freqs); 

if strcmpi(D.RecType,'acoustic') && isfield(D,'probename') && ~strcmpi(D.probename,'none'),
    ProbeGain = ProbeLoss(D.probename, 1e3*freqs);
    MG = MG - ProbeGain; % minus sign: we "undo" the probe transfer
end


%======= LOCAL =========
function [M, freqs] = local_NoisefromData(D, freqs)
% local_NoisefromData -return "evaluated" noisefloor as stored in ZWOAE data

if nargin < 2,
    freqs = (0:D.Nsam-1).'*D.df; % freq of all spectral components in kHz
end

%wrap frequencies to have them in interval 0...Nyquist
freqs = mod(freqs, D.fs);
freqs = min(freqs, D.fs-freqs);
M = polyval(D.Pnf, log(freqs)); %get the noise floor @ specified freqs
%NOTE: changing last line also requires changing fitNoise.m
