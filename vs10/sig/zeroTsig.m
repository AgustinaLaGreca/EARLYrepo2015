function S = zeroTsig(Fsam, Dur,t0);
% zeroTsig - zero-valued tsig.
%    Z=zeroTsig(Fsam, Dur,t0) generates a tsig object filled with zeros,
%    sampled at Fsam kHz, having a duration Dur ms and onset t0. 
%    Dur and t0 may be row vectors, in which case a multichannel tsig is 
%    generated.
%
%    See also tsig.

if nargin<3, t0=0; end

error(numericTest(Fsam,'scalar/rreal/positive', 'Sample frequency is '));
DT = 1/Fsam;
if (size(Dur,1)>1) || isempty(Dur),
    error('Dur argument must be row vector.');
end
if (size(t0,1)>1) || isempty(t0),
    error('t0 argument must be row vector.');
end
error(numericTest(Dur,'rreal/nonnegative', 'Dur argument is '));
error(numericTest(t0,'rreal', 't0 argument is '));

try,
    [Dur,t0] = SameSize(Dur,t0);
catch,
    error('Dur and t0 must have compitable channel counts.');
end
Nchan = size(Dur,2);
Nsam = round(Dur*Fsam);
W = cell(1,Nchan);
for ichan=1:Nchan,
    W{ichan} = zeros(Nsam(ichan),1);
end,
S = tsig(Fsam, W, t0);


