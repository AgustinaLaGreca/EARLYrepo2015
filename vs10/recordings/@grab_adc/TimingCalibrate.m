function TC=TimingCalibrate(G, flag);
% eventdata/TimingCalibrate - calibrate timing between analog output & input
%    TC=TimingCalibrate(G) performs a timing calibration of the Hardware and
%    returns the result in struct TC.
%
%    TimingCalibrate loads the appropriate SeqPlay-trig-2ADC circuit to
%    the D/A hardware (if needed) and calibrates the timing of its ADC 
%    recording re the analog output. The result is put in the TimingCalib
%    field of G in the form a a struct with fields:
%          Fsam: true sample frequency in Hz
%        Lag_ms: lag between D/A and A/D in ms
%    The circuit is only loaded when it hasn't been loaded already while
%    using the same sample rate.
%
%    TC=TimingCalibrate(G, 'force') always performs the calibration, even
%    if the same circuit is running at the same sample frequency since last
%    time.
%
%    The real work is delegated to DA_AD_TimingCalibrateRX6.
%
%    See also TimingCalibrateRX6, eventdata, grabevents, StimGUIAction,
%    SeqPlay.

if nargin<2, flag=''; end
persistent TClast IDlast


% check if the present "load" of the circuit has been calibrated before
CI = sys3CircuitInfo(G.Dev); % when sys3loadcircuit really loads a circuit, its gets a new ID
if ~isequal('RX6', G.Dev(1:3)),
    error(['Timing calibration not implemented yet for device ' G.Dev]);
end
if isequal(IDlast, CI.unique_id) && ~isequal('force',flag),
    TC = TClast;
    return;
end

TC = DA_AD_TimingCalibrateRX6(G.Dev, G.Hardware);
TClast = TC;
IDlast = CI.unique_id;








% 
