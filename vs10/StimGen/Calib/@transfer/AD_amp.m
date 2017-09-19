function A = AD_amp(T);
% transfer/AD_amp - waveform amplitude of A/D during transfer measurement
%    A = AD_amp(T) returns an array containing the peak amplitude of the
%    numbers obtained from the D/A converter during the transfer measurement.
%    (For the RX6 and RP2, this is the peak voltage of the FM sweep.)
%    The elements of A correspond to the frequencies in Frequency(T). 
%    
%    Note: the differences with Resp_level(T) are:
%      - Resp_level is in dB
%      - Resp_level referes to RMS value, not peak value
%      - a potentially non-unity value of the sensitivity factor of the 
%        response (see Transer/measure)
%
%    See Transfer, Transfer/Stim_level, Transfer/DA_amp.


if ~isfilled(T),
    error('Transfer object T is not filled.');
end

% first reconstruct the *stimulus* amplitude as used by Transfer/measure
L = Resp_level(T);
A = sqrt(2)*dB2A(L)*T.CalibParam.Sens.Resp; % dB RMS -> amp RMS -> Peak amp stim -> Peak amp A/D





