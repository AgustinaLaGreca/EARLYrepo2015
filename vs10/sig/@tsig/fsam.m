function Y = fsam(T, F)
% tsig/fsam - get/set sample rate of tsig object
%
%   fsam(T) or T.fsam is the sample rate [kHz] of T.
%
%   T=fsam(T,F) or T.fsam=F sets the sample rate to F kHz.
%   Note that changing the sample rate of T may affect the onset, too,
%   because the onset is always rounded toward the nearest integer number
%   of sample periods 1/F.
%
%   See tsig, "methodshelp tsig".

if nargin==1, % get
    Y = T.Fsam;
elseif nargin==2, % set
    T.Fsam = F;
    [T, Mess] = test(T); % let tsig constructor do the checking
    error(Mess);
    Y = T;
else, 
    error('Invalid #input arguments.');
end
    


