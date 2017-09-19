function X = ExactGate(X, Fsam, BurstDur, Tstart, RiseDur, FallDur, StartAmp, EndAmp);
% ExactGate - exactly timed gating using cos^2 ramps
%    X = ExactGate(X, Fsam, BurstDur, Tstart, RiseDur, FallDur) imposes 
%    cos^2 ramps onto waveform X. Inputs are
%          X: waveform. Must be single array.
%       Fsam: sample rate in Hz.
%   BurstDur: duration of burst in ms, measured from start-of-rise to 
%             end-of-fall. Sub-sample precision is respected.
%     Tstart: start of rising ramp in ms re first sample of X. Must be
%             nonnegative. Sub-sample precision is respected.
%    RiseDur: duration in ms of cos^2 rising ramp.
%    FallDur: duration in ms of cos^2 falling ramp.
%    StartAmp: Starting Amplitude of signal before gate (Optional; Default = 0) 
%    EndAmp: Amplitude of signal after X (Optional; Default = 0)
%
%    StartAmp & EndAmp can be used for creating gating signals if the
%    signal does not start/end from baseline but from another Burst with
%    Ampliture StartAmp/EndAmp. The window raises/gates the signal from StartAmp to the
%    amplitude of the second Burst max(X). The gating applied at the end of
%    X lowers/gates to a signal with amplitude EndAmp. 
%    Used for the STEP stimulus for
%
%    example. 
%    See also simplegate, toneStimSTEP.
if nargin < 7 
    StartAmp = [];
    EndAmp = [];
end

if ~isvector(X),
    error('X must be vector.');
end
[X, wasRow] = TempColumnize(X);

error(numericTest(Tstart, 'nonnegative/rreal/scalar', 'Tstart arg is'));
error(numericTest(RiseDur, 'nonnegative/rreal/scalar', 'RiseDur arg is'));
error(numericTest(BurstDur, 'nonnegative/rreal/scalar', 'BurstDur arg is'));
error(numericTest(FallDur, 'nonnegative/rreal/scalar', 'FallDur arg is'));

NsamTot = numel(X);
dt = 1e3/Fsam;
TotDur = dt*NsamTot;
if Tstart+RiseDur>TotDur,
    error('Rising ramp exceeds signal duration.');
end
if Tstart+BurstDur>TotDur+dt,
    Tstart, BurstDur, TotDur
    error('End of falling ramp exceeds signal duration.');
end
if Tstart+BurstDur-FallDur<0,
    error('Start of falling portion preceeds signal start.');
end

NsamHead = round((Tstart+RiseDur)/dt); % # samples in the heading silence plus rising portion
t_re_tstart = dt*(0:NsamHead-1)-Tstart;
if RiseDur>0,
    if ~(isempty(StartAmp) || isempty(EndAmp))
        A1 = StartAmp;
        A2 = max(X);
        a = A1/A2;
        b = 1-a;
        Wrise = a+b*sin(2*pi*t_re_tstart/(4*RiseDur)).^2.'; % generate and 
        X(1:NsamHead) = Wrise.*X(1:NsamHead);  % ... apply rising window
    else
        Wrise = sin(2*pi*t_re_tstart/(4*RiseDur)).^2.'; % generate and 
        X(1:NsamHead) = Wrise.*X(1:NsamHead);  % ... apply rising window
    end
end
X(find(t_re_tstart<0)) = 0; % set pre-rise samples to zero.

Tfall = Tstart+BurstDur-FallDur; % start of falling portion
isamTail = (1+ceil(Tfall/dt):NsamTot).'; % sample indices of falling portion
t_re_tfallstart = (isamTail-1)*dt - Tfall; % time re start of fall
if FallDur>0,
    if ~(isempty(StartAmp) || isempty(EndAmp))
        A2 = max(X);
        A3 = EndAmp;
        a = A3/A2;
        b = 1-a;
        Wfall = a+b*cos(2*pi*t_re_tfallstart/(4*FallDur)).^2; % generate and ...
        X(isamTail) = Wfall.*X(isamTail); % apply fall window
    else
        Wfall = cos(2*pi*t_re_tfallstart/(4*FallDur)).^2; % generate and ...
        X(isamTail) = Wfall.*X(isamTail); % apply fall window
    end
end
Tend = Tstart+BurstDur; % end of falling portion
isamPost = (1+ceil(Tend/dt):NsamTot); % sample indices of post-fall era: ...
X(isamPost) = 0; % ... zero them


% restore original 'orientation'
if wasRow, X = X.'; end









