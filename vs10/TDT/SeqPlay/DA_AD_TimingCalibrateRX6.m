function TC=DA_AD_TimingCalibrateRX6(Dev, Hardware, doPlot);
% DA_AD_TimingCalibrateRX6 - calibrate timing between analog input & outputof RX6   
%    S=DA_AD_TimingCalibrateRX6(Dev, Hardware)
%    measures the true delay between AD and DA of the RX6 by playing a 
%    noise burst and evaluating the crosscorrelation between played and 
%    recorded waveforms. It is assumed that the appropriate SeqPlay
%    circuit has been loaded using the correct sample rate.
%
%    The result of the calibration is returned in struct S whose fields are
%          Fsam: true sample frequency in kHz
%        Lag_ms: lag in ms between D/A and A/D
%      PeakCorr: peak normalized correlation between played & recorded
%                waveforms. Should be >0.9. Lower values indicate wiring
%                problems.
%
%    The default flaf is Default flag is ''. That is, 
%    TimingCalibrateRX6(Fsam, Postfix, '') is the same as 
%    TimingCalibrateRX6(Fsam, Postfix).
%
%    uses bit CalStampBit for acquiring the pulse train and afterwards sets 
%    the timestamp bit to SpikeStampBit for acquiring spike times.
%    Defaults are CalStampBit=3; SpikeStampBit=0.
%
%    See also TimingCalibrate, SeqPlay, sys3loadCircuit.

if nargin<3, doPlot=0; end
Fsam = sys3CircuitInfo(Dev, 'Fsam')*1e3; % sample rate in Hz
Nsam = 200; Ups = 500;
Dur = 1e3*Nsam/Fsam; % dur of calib signal in ms

% suppress sound delivery
% SetAttenuators(experiment,'max','max', Hardware);
SetAttenuators(experiment,'max',nan, Hardware); %cahgned by EVE 13/08/2015

% set duration of relay activation
sys3setpar(max(200, 1.5*Dur), 'CalibDur_ADC', Dev);

% prepare & upload wave
Wav = 0.5*randn(Nsam,1);
Wav = [Wav(:); zeros(Nsam,1)];
seqplayupload({Wav},{});
seqplaylist(1,1);
% enable event grabbing by the Seqplay circuit
seqplaygrab_adc on;
% reset the ADC buffers
feval(Hardware.ResetADC_trigger); % empty the ADC buffers
feval(Hardware.CalibADC_trigger); % activate the relay connecting DAC to ADC
pause(0.05); % small delay for the relay to kick in
% play & record & upsample 
seqplaygo;
seqplayDAwait;
Rec = sys3read('ADC_1', 2*Nsam, Dev);
Wav = resample(Wav,Ups,1);
Rec = resample(Rec,Ups,1);

% 24/03/2014: TVP problem with DC offset of -.25V on input channel 1 (internal
% problem?) --> to be checked
% remove DC
Wav = Wav - mean(Wav);
Rec = Rec - mean(Rec);

% crosscorrelate upsampled Wav & Rec and estimate lag by locating maximum
dt = 1e3/Fsam/Ups; % sample period in ms of upsampled waveforms
[XX,lag]=xcorr(Rec,Wav, 'coeff');
lag = lag*dt;
[PeakCorr imax] = max(XX);
if doPlot || (PeakCorr<0.85),
    f8; subplot(1,2,1); dplot(dt, Wav); xdplot(dt, Rec,'r');
    puttext('northwest', 'A/D problems??', 'fontsize', 18);
    subplot(1,2,2); dplot(dt*[1 min(lag)], XX); 
    error('DA->AD timing calibration of RX6 failed.');
end
ir = imax+(-10:10);
lag = lag(ir);
XX = XX(ir);
%subplot(2,1,1); plot(lag ,XX,'.-');
dXX = diff(XX);
lag = (lag(1:end-1)+lag(2:end))/2;
%subplot(2,1,2); plot(lag,dXX); grid on
Lag_ms = interp1(dXX,lag,0,'linear');
% wait till calib relay has flipped back
while 1,
    if ~sys3getpar('RelayStatus_ADC', Dev), 
        break;
    end
end
TC = CollectInStruct(Fsam, Lag_ms, PeakCorr);


%=====obsolete=====
% % we don't have the right relays to switch off the regular inputs to the 
% % analog IN (e.g. neural signal, microphone) during calibration, so we'll
% % just use the DAC-versus-digital delay of the ()
% if ~isequal('RX6', G.Dev),
%     error('Timing calibration of non-RX6 device not implemented');
% end
% TC = TimingCalibrateRX6(G.Fsam/1e3);
