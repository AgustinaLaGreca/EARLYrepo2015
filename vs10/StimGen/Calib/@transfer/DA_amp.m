function A = DA_amp(T);
% transfer/DA_amp - amplitude of waveform sent to D/A transfer measurement
%    A = DA_amp(T) returns an array containing the peak amplitude of the
%    numbers sent to the D/A converter during the transfer measurement.
%    (For the RX6 and RP2, this is the peak voltage of the FM sweep.)
%    The elements of A correspond to the frequencies in Frequency(T). 
%    
%    Note: the differences with Stim_level(T) are:
%      - Stim_level is in dB
%      - Stim_level referes to RMS value, not peak value
%      - a potentially non-unity value of the sensitivity factor of the 
%        stimulus (see Transer/measure)
%
%    See Transfer, Transfer/Stim_level, Transfer/AD_amp.


if ~isfilled(T),
    error('Transfer object T is not filled.');
end

% first reconstruct the *stimulus* amplitude as used by Transfer/measure
L = Stim_level(T);
A = sqrt(2)*dB2A(L)*T.CalibParam.Sens.Stim; % dB RMS -> amp RMS -> Peak amp stim -> Peak amp D/A





