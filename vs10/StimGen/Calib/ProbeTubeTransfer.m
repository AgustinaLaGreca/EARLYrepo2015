function [Freq, Magn, Phase]=ProbeTubeTransfer(TRFprobe, TRFtm, ProbeSens, TmSense, FN);
% ProbeTubeTransfer - probe transfer function from microphone data
%   [Freq, Magn, Phase] = ProbeTubeTransfer(TRFprobe, TRFtm, ProbeSens, TmSense, FN)
%   Inputs:
%       TRFprobe, TRFtm: either the structs as returned by MeasureTransfer
%                        or names of mat files holding thses structs.
%    ProbeSens, TmSense: sensitivity (V/Pa) of probe and tm microphone,
%                        respectively.
%                    FN: optional filename to save this transfer function. 
%   Outputs:
%       Freq: frequency in Hz
%       Magn: tm-to-probe gain in dB (loss=negative)
%      Phase: tm-to-probe phase increment in cycle (lag=negative)
if nargin<5, FN=''; end % default: don't save

Nfreq = 5e4; % # samples

TRFprobe = local_read(TRFprobe);
TRFtm = local_read(TRFtm);
if ~isequal(TRFprobe.Freq, TRFtm.Freq),
    error('Transfer functions of probe and TM have incompatible frequency axes.');
end
% trfplot({TRFtm, TRFprobe}); 

% ---compute TM->Probe TRF & iinclude net lag in frequency domain
comp_lag = TRFprobe.Lag_ms - TRFtm.Lag_ms;
comp_phase = 1e-3*TRFprobe.Freq.*comp_lag;
phasor = exp(-2*pi*i*comp_phase);
gainor = TmSense/ProbeSens;
TRF = gainor.*phasor.*TRFprobe.TRF./TRFtm.TRF;
Magn = A2dB(abs(TRF));
Phase = cunwrap(cangle(TRF));
Freq = TRFprobe.Freq;
%$ decimate or interpolate
fr = logispace(min(Freq),max(Freq),Nfreq);
Magn = interp1(Freq, Magn, fr);
Phase = interp1(Freq, Phase, fr);
Freq = fr;
local_save(Freq, Magn, Phase, FN);
% f2; trfplot(TRFprobe)


%===========================================
function TRF = local_read(TRF);
% if TRF is filename, read the file
if ~ischar(TRF), return; end
DD = fullfile(ZWOAEdatadir, 'calibration');
TRF = fullfilename(TRF, DD, '.trf');
if ~exist(TRF,'file'),
    error(['File ''' TRF ''' not found.']);
end
qq=load(TRF);
FN = fieldnames(qq);
if numel(FN)>1,
    error(['File ''' TRF ''' contains multiple variables, not just one transfer function.']);
end
TRF = qq.(FN{1});

function local_save(Freq, Magn, Phase, FN);
if isempty(FN), return; end
TRF = CollectInStruct(Freq, Magn, Phase);
DD = fullfile(ZWOAEdatadir, 'calibration');
FN = fullfilename(FN, DD, '.probeloss');
save(FN, 'TRF', '-mat');











