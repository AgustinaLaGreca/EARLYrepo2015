function S=TimingCalibrateRX6(Fsam, Postfix, flag, CalStampBit, SpikeStampBit);
% TimingCalibrateRX6 - calibrate timing between analog output & digital input of RX6   
%    S=TimingCalibrateRX6(Fsam) loads the RX6SeqPlay-triggering circuit to the RX6 (if needed)
%    and calibrates its time stamping. The result is returned in struct S 
%    whose fields are
%       Fsam: true sample frequency in kHz
%        Lag: lag between D/A and digital input in 
%    Details: detailed parameters of the calibration.
%    The circuit is only loaded when it hasn't been loaded before using the
%    same sample rate.
%
%    S=calibrateRX6(Fsam, '-foo') loads the circuit RX6SeqPlay-foo instead 
%    of the default circuit RX6SeqPlay-trig-2ADC. 
%
%    TimingCalibrateRX6(Fsam, Postfix, 'force') forces loading of the circuit.
%    The default flaf is Default flag is ''. That is, 
%    TimingCalibrateRX6(Fsam, Postfix, '') is the same as 
%    TimingCalibrateRX6(Fsam, Postfix).
%
%    S=TimingCalibrateRX6(Fsam, Postfix, flag, CalStampBit, SpikeStampBit, CalibTrigger)
%    uses bit CalStampBit for acquiring the pulse train and afterwards sets 
%    the timestamp bit to SpikeStampBit for acquiring spike times.
%    Defaults are CalStampBit=3; SpikeStampBit=0.
%
%    See also TimingCalibrate, SeqPlay, sys3loadCircuit.

persistent last_S

if nargin<2, Postfix = '-trig-2ADC'; end
if nargin<3, flag = ''; end
if nargin<4, CalStampBit=3; end
if nargin<5, SpikeStampBit=0; end

RecycleCircuit = ~isequal('force', flag); % default is to reuse loaded circuit

if Fsam>131,
    error('Sample rate may not exceed 130.21 kHz; circuit will fail.');
end

% stimulus parameters
Vbase = -1.5;
Atrain = 3; % V pulse amp
Ftrain = 127; % Hz pulse rate of train
Npulse = 15; % # pulses

% Initialize sequenced play; generate pulse waveform & upload to RX6
Fsam = RX6sampleRate(Fsam);
[SPinfo, Recycled] = seqplayinit(Fsam, 'RX6', Postfix, RecycleCircuit);
sys3trig(3,'RX6'); % reset event buffer
if Recycled && ~isempty(last_S) && isequal(Fsam, last_S.Fsam),
    % same circuit, running @ same sample freq, need no retesting
    S = last_S;
    return;
end
% enable event grabbing
seqplaygrab(1);
% generate pulse train
Fsam = SPinfo.Fsam; % exact (TDT-compliant) sample rate in kHz
dt = 1/Fsam; % ms sample period
Nper = round(1e3*Fsam/Ftrain); % #samples in single pulse interval
dt_train = dt*Nper; % ms inter-pulse interval
Nramp = 5; % duration of ramp expressed
Prepulse = linspace(0,Vbase,Nramp*Nper).'; % gradual ramp toward negative prepotential
Pulse = [Atrain*[1]; Vbase*ones(Nper-1,1)]; % pulse waveform: single-sample positive excursion from neg prepotential
PostPulse = flipud(Prepulse); % gradual ramp from negative prepotential to zero
seqplayupload({Prepulse Pulse PostPulse},{0*Prepulse 0*Pulse 0*PostPulse});
seqplaylist([1 2 3], [1 Npulse 1],[1 2 3], [1 Npulse 1]);
%dplot(dt,Pulse); pause; delete(gcf)

% set parameters of trigger circuit
StimDur = (Npulse+2*Nramp)*dt_train; % ms total stimulus (ramps+pulse train) duration
sys3setpar(1.5*StimDur,'CalibDur',SPinfo.Dev); % time in ms that DAC-out is connected to Digital-in-1
sys3setpar(CalStampBit,'TstampBit',SPinfo.Dev); % time stamping from bit 1 (patch-panel connector A2)

% activate the circuit: 1) temp connect DA-> bit 1; 2) play train
sys3trig(2, SPinfo.Dev); % triggers connection DAC-out <-> Digital-in-1
pause(0.1); % circuit needs some time to realize switching
seqplaygo; % start D/A
seqplayDAwait; % wait for D/A conversion to finish
pause(0.1)
% read timing info from buffers & convert to event times (ms). 
EventTick = sys3read('EventTick', 'EventCount', SPinfo.Dev, 0, 'I32');
EventTstamp = 1e3*sys3read('EventTstamp', 'EventCount', SPinfo.Dev, 0); % ms
Tevent = dt*EventTick + EventTstamp;
% Weed out any false triggers
Evdif = diff(Tevent); % intervals between consecutive events
iwrong = find((Evdif<0.1*dt_train) | (Evdif>1.8*dt_train)); % events following one another too fast or slowly
Tevent(1+iwrong)=[];
Nstamp = length(Tevent);

% To protect the digital input, wait until the relay that connects analog ..
while sys3getpar('RelayStatus', SPinfo.Dev), end % ..output & digital input is disabled
% Change the digital input for time stamping of the real spikes 
sys3setpar(SpikeStampBit,'TstampBit',SPinfo.Dev);

if ~isequal(Nstamp,Npulse),
    plot(EventTick, EventTick*0, '.');
    error('#pulses recorded ~= #pulses played. Check wiring!');
end

Tplay = Nramp*dt_train + dt_train*(0:Npulse-1); % "true" time of pulses when played
ir =3:Npulse;
P = polyfit(Tplay(ir),Tevent(ir),1);
Lag = P(2);
dSlope = P(1)-1;
NsamLag = Lag/dt;
LagFrac = rem(NsamLag, 1);
LagFracMicrosec = dt*1e3*LagFrac;

% plot(Tplay,Tevent,'.'); ff
% plot(diff(Tevent)-dt_train, '*'); ff

Lag_us = 1e3*Lag; % lag in us
Details = CollectInStruct(Tplay, Tevent, Fsam, dt, P, dSlope, Lag, Lag_us, NsamLag, LagFrac,LagFracMicrosec);
S = CollectInStruct(Fsam, Lag_us, Details);
% store for future use
last_S = S;
seqplaygrab(0); % disable event grabbing
sys3trig(3,'RX6'); % reset event buffer for measurement





