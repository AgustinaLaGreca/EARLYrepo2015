function C = NewMeasureTransfer(Dev, ichan_out, ichan_in, Frange, Fsam, FMspeed, Amp);
% NewMeasureTransfer - measure transfer function using RP2 DAC & ADC
%   C = measureTransfer(Dev, ichan_out, ichan_in, Frange, Fsam, FMspeed, Amp);
%   Fsam in kHz (will be "rounded" to nearest TDT sample rate)
%   Flow, Fhigh in Hz (!). Defaults 100 Hz and 0.45 times sample rate.
%   Amp in Volt; default 0.1 V. Amp may also be function handle Amp(Freq).
%   FMspeed in octaves/second; default: 1 oct/s
%   ichan_out: DAC channel: 1 [default] or 2.
%   ichan_in: ADC channel: 1 [default] or 2.
%
%   measureTransfer uses an FMsweep.
%
%   See also calibstim.

WinWidth = 0.05; % octave-wide smoothing

if nargin<1, Dev = []; end; % default value; see below
if nargin<2, ichan_out=1; end
if nargin<3, ichan_in=1; end

if nargin<4, Frange = []; end % default value; see below
if nargin<5, Fsam = []; end % Volt default stim ampl
if nargin<6, FMspeed=1; end % default 1 oct/s
if nargin<7, Amp=1; end % default 0.1 Volt

if isempty(Fsam), Fsam = 50; end % 50 kHz default
Fsam = TDTsampleRate(Dev, Fsam); % exact sample rate in kHz for this device

if isempty(Frange), Frange = [50 400*Fsam]; end % from 50 Hz to 40% of the sample freq
[Flow, Fhigh] = DealElements(Frange); % extreme freqs in Hz
% smoothing window
NsamOctWin = 1+2*round(1e3*Fsam*WinWidth*FMspeed/2); % # stimulus samples spanning WinWidth octaves during FM sweep; ...
OctWin = hanning(NsamOctWin); % ... make sure its an odd number
OctWin = OctWin/sum(OctWin); % normalize

% load circuit & select channnels
Fsam = sys3loadCircuit('AcoustCalibRP2', Dev, Fsam);
sys3setpar(ichan_out, 'DACchan', Dev);
sys3setpar(ichan_in, 'ADCchan', Dev);

% upload waveform
NsamBuf = sys3ParTag(Dev, 'PlayBuf', 'TagSize');
% does the whole sweep fit in the play buffer? If not, divide in chunks.
S=calibstim(1e3*Fsam, Flow, Fhigh, FMspeed, Amp);
% do noise calib to get estimate of lag
Rec=local_play_rec(Dev, S, 0, 'N');
[MaxNoiseCor, NsamLag] = maxcorr(Rec, S.N); % sample shift  yielding highest corr
Lag_ms = NsamLag*S.dt;

% now that we now the lag, we can compute # chunks
Nchunk = ceil((S.Nsam+NsamLag)/NsamBuf);
EdgeFreq = logispace(Flow, Fhigh, Nchunk+1);

TRF = [];
Freq = [];
% perform calib chunk by chunk
for ichunk =1:Nchunk,
    freq0 = EdgeFreq(ichunk); freq1 = EdgeFreq(ichunk+1);
    % generate & upload waveform
    S=calibstim(1e3*Fsam, freq0, freq1, FMspeed, Amp);
    Rec=local_play_rec(Dev, S, NsamLag);
    Rec=localRamp(Rec(:),S.NsamRamp); % avoid boundary artifacts by ramping uninformative ends
    % heterodyne the recorded response to zero-freq by compensating for the stimulus freq
    Rec = fft(Rec);
    Rec = local_limitSpec(Rec, 1e3*Fsam, 0.7*freq0, 1.3*freq1);
    Rec = ifft(Rec);
    Rec = Rec.*conj(S.Sweep)/Amp.^2*0.5;
    Rec = localRamp(Rec(:),S.NsamRamp); % avoid boundary artifacts by ramping uninformative ends
    isweep = S.NsamRamp+2:S.Nsam-S.NsamRamp; % indices of proper sweep w/o ramps
    TRF = [TRF; Rec(isweep)];
    Freq = [Freq; S.InstFreq(isweep)];
end
%dsize(TRF, OctWin);
TRF = localWin(TRF,OctWin); % remove fast ripples in transfer function 
% figure;
% subplot(2,1,1);
% plot(Freq/1e3, a2db(abs(TRF))); 
% xlog125([0.1 Fsam/2]);
% subplot(2,1,2);
% plot(Freq/1e3, unwrap(angle(TRF))/2/pi); 
% xlog125([0.1 Fsam/2]);

C = CollectInStruct(Fsam, Amp, FMspeed, ichan_out, ichan_in, NsamOctWin, EdgeFreq, NsamLag, Lag_ms, MaxNoiseCor, Freq, TRF);

%======================
function Rec=local_play_rec(Dev, S, NsamLag, FN);
if nargin<3, NsamLag = 0; end % default: no lag correction
if nargin<4, FN = 'Sweep'; end % default: play sweep
W = real(S.(FN)); Nsam = numel(W);
sys3write([W; zeros(NsamLag,1)], 'PlayBuf', Dev);
sys3trig(1,Dev); % reset buffer indices of both play buffer & record buffer
sys3run(Dev);
sys3waitfor('NsamStored', '>=', Nsam+NsamLag, Dev);
Rec = sys3read('RecordBuf', Nsam, Dev, NsamLag); % recording is retarded by NsamLag samples re stimulation
sys3halt(Dev);

function W=localRamp(W,NsamRamp);
% supply cos^2 onset & offset ramps
Ramp = sin(linspace(0,pi/2,NsamRamp).');
W(1:NsamRamp) = Ramp.*W(1:NsamRamp);
W(end-NsamRamp+(1:NsamRamp)) = flipud(Ramp).*W(end-NsamRamp+(1:NsamRamp));

function RS = local_limitSpec(S, FS, F0, F1);
% restrict complex spectrum S to positive freqs between F0 and F1
df = FS/numel(S); % freq spacing
i0 = 1+round(F0/df); i1 = 1+round(F1/df); % freq -> spectral index
RS = 0*S;
RS(i0:i1) = 2*S(i0:i1); % factor 0.5 to compensate for loss of neg freqs


function TRF = localWin(TRF,OctWin);
NdelayOctWin = round((numel(OctWin)-1)/2); % lag introduced by convoving with OctWin 
reps = ones(NdelayOctWin,1);
TRF = [reps*TRF(1); TRF; reps*TRF(end)]; % repeat first & last sample
TRF = conv(TRF, OctWin); % apply "lowpass" window, throwing out rapid flutuations
Ncut = 2*NdelayOctWin; % samples to cut from conv
TRF = TRF(Ncut+1:end-Ncut); % remove tails introduced by repeating first&last samples, and by convolution itself




