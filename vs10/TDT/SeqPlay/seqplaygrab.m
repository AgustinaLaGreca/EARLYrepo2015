function y = SeqplayGrab(G);
% SeqplayGrab - set or get grabbing status of seqplay circuit
%   SeqplayGrab returns the grabstatus of the currently loaded seqplay
%   circuit. When loading the circuit, the default status is 0, which
%   means that event grabbing is disabled. 
%
%   SeqplayGrab(1) or SeqplayGrab('on') enables event grabbing.
%
%   SeqplayGrab(0) or SeqplayGrab('off') enables event grabbing.
%
%   See also SeqPlay, Seqplayinit, SeqplayInfo.


P = private_seqPlayInfo; % info shared by seqplayXXX

if nargin<1, % get
    y = sys3getpar('DoGrab', P.Dev);
else, % set
    if isequal('on', G), G = 1; 
    elseif  isequal('off', G), G = 0; 
    end
    sys3setpar(G, 'DoGrab', P.Dev);
end


