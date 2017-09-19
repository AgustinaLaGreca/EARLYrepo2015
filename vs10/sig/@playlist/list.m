function [Y,Z] = list(P,L)
% playlist/list - get/set the list of playlist object
%   list(P) or P.list returns the waveform indices and repetition counts of
%   playlist object P in the format described below.
%   [iWave, Nrep] = list(P) returns the waveform indices and rep counts in
%   separate row vectors.
%
%   P=list(P,[iWav; Nrep]) or P.list=[iWav; Nrep] specifies the
%   waveform indices and repetition counts of listplay object P.
%   iWave and Nrep are row vectors containing the waveform indices and rep
%   counts, respectively.
%
%   See playlist, playlist/iWave, playlist/Nrep.

if nargin==1, % get
    Y = P.iplay;
    Z = P.Nrep;
    if nargout<2, % put in single matrix
        Y = [Y;Z];
    end
elseif nargin==2, % set
    if isvoid(P),
        error('Cannot set list of void playlist object.')
    else, % split matrix into iWave & Nrep vectors (see help), test, & set.
        if ~isequal(2,size(L,1)), 
            error('List spec must be 2xN integer matrix.'); 
        end
        iWave = L(1,:); Nrep = L(2,:);
        error(numericTest(iWave,'posint/noninf/real', 'iWave spec is '));
        if any(iWave>length(P.Waveform)),
            error('Waveform index in list exceeds # waveforms stored.');
        end
        error(numericTest(Nrep,'nonnegint/noninf/real', 'Nrep spec is '));
        P.iplay = iWave;
        P.Nrep = Nrep;
        Y = P;
    end
else,
    error('Wrong # input args')
end









