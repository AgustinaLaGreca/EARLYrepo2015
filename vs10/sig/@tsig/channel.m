function T = channel(T, k, S)
% tsig/channel - select or replace channels of tsig object
%
%   channel(T,K) or T(K) is T restricted to channels K,. K is a row vector
%   of positive integers <= nchan(T).
%
%   T=channel(T,k,S) or T(k)=S replaces the kth channel of T by S.
%   k must be scalar index; S must be single-channel. 
%
%   T=channel(T,K,[]) or T(k)=[] removes channel(s) K from S.
%
%
%   See tsig/nchan, tsig/subsasgn, tsig/horzcat.

if nargin<2, return; end

Nc = length(T.Waveform);

if nargin<3, % get
    if isequal(':', k), k = 1:Nc; end
    error(numericTest(k,'posint/noninf','Channel index is '));
    if ~isvector(k) || size(k,1)>1,
        error('Index K must be scalar or row vector.');
    end
    if any(k>Nc),
        error('Channel index exceeds number of channels in tsig.');
    end
    T.Waveform = T.Waveform(k);
    T.t0 = T.t0(k);
else, % set
    if isequal([],S), % remove channels
        if isempty(k), return; end; % T([])=S allowed, just as in ordinary Matlab assigments
        error(numericTest(k,'rreal/posint','Channel index k in subscripted assignment T(k)=S is '));
        if size(k,1)>1, error('Channel index k in assignment  T(k)=[]  must a row vector.'); end
        if k>nchan(T), error('Index k in   T(k)=[]  exceeds number of channels of T.'); end
        T.Waveform(k) = [];
        T.t0(k) = [];
        if isempty(T.t0),
            error('It is not allowed to remove all channels from a tsig object as in T(:)=[].');
        end
    elseif isTsig(S),
        if isempty(k), return; end; % T([])=S allowed, just as in ordinary Matlab assigments
        error(numericTest(k,'rreal/posint','Channel index k in assignment T(k)=S is '));
        if size(k,1)>1, error('Channel index k in assignment  T(k)=S  must a row vector.'); end
        if ~isequal(1,nchan(S)) && ~isequal(length(k),nchan(S)),
            error('Incompatible number of channels of LHS and RHS in assignment  T(k)=S.');
        end
        if ~isequal(T.Fsam,S.Fsam),
            error('In an tsig assignment T(k)=S,  S and T must have the same sample frequency.');
        end
        T.Waveform(k) = S.Waveform;
        T.t0(k) = S.t0;
    else,
        error('In an tsig assignment T(k)=S,  S must be [] or a tsig object.');
    end
end
[T, Mess] = test(T); 
error(Mess);






