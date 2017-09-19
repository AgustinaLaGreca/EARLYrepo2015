function [SP, dt, tAP] = JLspont(S, Fh, doReject625Hz);
% JLspont - spontaneous activity preceding stimulation of JLbeat data
%    [SP, dt] = JLspont(S) returns column array SP containing the
%    spontaneous activity of condition S. S may be JLbeat output or
%    JLbeatStats output. For arrays S, SP will be a matrix whose colmns are
%    the respective spontaneous waveforms. dt is sample period in ms.
%    Waveforms are highpass filtered from 20 Hz.
%    
%    JLspont(S, Fh) uses highpass filter with cutoff Fh (Hz).
%
%    JLspont(S, Fh, 1) rejects 8-Hz-wide band around 625 Hz
%
%    See also JLspontAna

if nargin<2, Fh=[]; end
if nargin<3, doReject625Hz=0; end

if isempty(Fh), Fh = 20; end % 20 Hz default

if isstruct(S) && numel(S)>1,
    SP = [];
    for ii=1:numel(S),
        [sp, dt, tAP{ii}] = JLspont(S(ii), Fh);
        SP = [SP, sp];
    end
    return;
end

% --------single S from here---------
CFN = mfilename;
CPM = {S.UniqueRecordingIndex, Fh, doReject625Hz};
YY = getcache(CFN, CPM);
if ~isempty(YY),
    SP = YY.SP;
    dt = YY.dt;
    tAP = YY.tAP;
    return;
end
% compute
StimOnset = 500; % ms start of beat stim
D = readTKABF(S);
dt = D.dt_ms;
R0 = D.AD(1).samples; clear D;
NsamSpont = 1+round(StimOnset/dt);
[dum, tAP] = APtruncate2(R0(1:NsamSpont), S.APthrSlope, dt, [S.APwindow_start S.APwindow_end]);
R0 = APtruncate2(R0, S.CutoutThrSlope, dt, [S.APwindow_start S.APwindow_end]);
[B,A] = butter(5, 2e-3*Fh*dt ,'high');
R = filter(B,A,R0);
%dplot(dt, R); xdplot(dt, R0, 'r'); pause; delete(gcf);
SP = R(1:NsamSpont);
SP = SP(:);
if doReject625Hz,
    SP = local_reject_625(SP,dt, rem(S.UniqueRecordingIndex,61^5));
end
% store in cache
putcache(CFN, 4e3, CPM, CollectInStruct(SP,dt, tAP));

%=============
function R = local_reject_625(R,dt,Rseed);
Nsam = numel(R);
df = 1/(Nsam*dt); % freq spacing in kHz
freq = Xaxis(R, df);
Rsp = fft(R);
Rsp(round(Nsam/2):end) = 0; % discard "negative freqs"
Rsp(1) = Rsp(1)/2;
Rsp = 2*Rsp; % anticipate iift; compensate for discarding the neg freqs
ilo = find(betwixt(freq, 0.601, 0.617));
istop = find(betwixt(freq, 0.617, 0.631));
ihi = find(betwixt(freq, 0.633, 0.647));
MNP = mean(abs(Rsp([ilo;ihi])).^2); % mean power in flanking bands
STP = mean(abs(Rsp([ilo;ihi])).^2); % power std in flanking bands
SetRandState(Rseed);
NewStopAmplitudes = sqrt(abs(MNP + STP*randn(size(istop)))); % subsitute amplitudes
Rsp(istop) = NewStopAmplitudes.*Rsp(istop)./abs(Rsp(istop));
%dplot(df, p2db(abs(Rsp))); pause; 
R = real(ifft(Rsp));










