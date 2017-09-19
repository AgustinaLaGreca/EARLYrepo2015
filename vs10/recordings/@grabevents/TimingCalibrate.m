function TC=TimingCalibrate(G);
% grabevents/TimingCalibrate - calibrate timing between analog output & digital input
%    TC=TimingCalibrate(G) performs a timing calibration of the Hardware and
%    returns the result in struct TC.
%
%    TimingCalibrate loads the appropriate SeqPlay-triggering circuit to
%    the D/A hardware (if needed) and calibrates its time stamping re the
%    analog output. The result is put in the TimingCalib field of ED in
%    the form a a struct with fields:
%       Fsam: true sample frequency in Hz
%        Lag: lag between D/A and digital input in
%    Details: detailed parameters of the calibration.
%    The circuit is only loaded when it hasn't been loaded already while
%    using the same sample rate.
%
%    These fields are used to obtain the information about the stimulus and
%    the hardware needed for the timing calibration. The timing calibration
%    itself is delegated to functions like TimingCalibrateRX6, which play a
%    pulse train over the D/A and record it over the digital input.
%    Timing calibration is only performed when a new circuit and/or sample
%    rate is needed than last time; otherwise the outcome of the previous
%    test is used.
%
%    See also TimingCalibrateRX6, eventdata, grabevents, StimGUIAction,
%    SeqPlay.

% suppress sound delivery
SetAttenuators(experiment(),'max','max', G.Hardware);

% delegate to hardware-specific function
TCfun = fhandle([mfilename G.Hardware.DAC]); % @timingCalibrateRX6 etc

okay=0;
for ii=1:3, % calibration sometimes fails for no good reason. Give it 3 shots before giving up
    try,
        TC = feval(TCfun, G.Fsam/1e3, G.Hardware.SeqplayPostFix,'force', ...
            G.Hardware.CalibTimeStampBit, G.Hardware.SpikeTimeStampBit);
        okay = 1;
    catch,
        warning(lasterr);
    end
end
if ~okay, rethrow(lasterror); end







