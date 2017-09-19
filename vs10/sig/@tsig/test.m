function [T, mess] = test(T);
% tsig/test - test contents of tsig object.
%
%   [T, Mess]=test(T) tests whether the content of tsig object T is valid 
%   by first extracting the contents of T, and then passing them to the 
%   tsig constructor. Any error message from tsig are caught and returned 
%   in Mess.
%
%   See tsig.

mess = '';
try, T = tsig(T.Fsam, T.Waveform, T.t0, T.Internal);
catch, mess = lasterr;
end










