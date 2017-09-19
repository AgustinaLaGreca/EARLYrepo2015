function m = MaxAbsDA(T);
% transfer/MaxAbsDA - max absolute value that D/A converter can produce
%    MaxAbsDA(T) returns the max abs value of the samples that can be sent
%    to the D/A converter that was used to measure T.
%
%    See Transfer, Transfer/measure


if ~isfilled(T),
    error('Transfer object T is not filled.');
end
CP = [T.CalibParam]; % handle arrays
Amp = [CP.Amp];
m = [Amp.MaxAbsDA];



