function Freq= frequency(T);
% transfer/frequency - frequency of response during transfer measurement
%    frequency(T) returns an array containing the frequencies in Hz as used
%    during the measurement of T. 
%
%    See Transfer, Transfer/measure, Transfer/resp_amp.


if ~isfilled(T),
    error('Transfer object T is not filled.');
end

Freq = T.Freq;




