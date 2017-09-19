function [Dt, binDt]=CalibDelay(EXP, Fsam, DAchan, Nchan);
% calibDelay - effective latencies of DA conversion
%   [Dt, binDt]=calibrate(EXP, Fsam, chanStr) returns the latencies of
%   D/A conversion computed from the calibration data. Outputs are
%         Dt: global group delay [ms] of the calibration data. Estimate of 
%             the electrical+acoustical delay. If a wideband waveform is 
%             converted at time zero, the mean arrival time of its 
%             acoustical energy is at Dt. This group delay should *not* be
%             compensated when computing stimulus waveforms (see Calibrate) 
%             The Dt parameter, however, is important at analysis time, 
%             when interpreting spike data and their latencies and group 
%             delays re the acoustical waveform.
%  binCompDt: zero if this DAchannel is lagging the other in terms of Dt,
%             minus the difference if this DAchannel is the leading one.
%             Delaying the waveform of this DAchannel by binDt, will 
%             compensate the interaural asymmetry in latency caused by 
%             acoustical or electrical differences in the D/A-to-eardrum 
%             path. When only one channel is active (monaural calibration 
%             data), binCompDt is zero by convention.
% Nchan: Number of channels used for this stimulus (monaural: 1 Binaural: 2)  
%
%   Note that both Dt and binDt may vary with Fsam. See Calibrate for
%   details.
%          
%   See also Calibrate.

% delegate to calibrate
[dum, dum, Dt, binDt]=calibrate(EXP, Fsam, DAchan,[]);
if Nchan == 1
    binDt = 0;
end









