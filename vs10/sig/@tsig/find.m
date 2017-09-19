function t = find(T,K,varargin)
% tsig/find - find time instants of nonzero values.
%   TC=find(T) returns a cell array TC with TC{ichan} a column vector 
%   containing the time [ms] instants of the nonzero entries of the 
%   waveform in T(ichan).
%
%   TC=find(T,K) returns at most K time values per channel. K may be a 
%   T.nchan-long row vector, in which case at most K(ichan) are returned 
%   for channel ichan.
%
%   TC=find(T,K, 'first') is the same as TC=find(T,K) 
%   
%   TC=find(T,K, 'last') returns at most the K last time values
%   corresponding to the nonzero entries of T. Again, K is applied in
%   a channelwise fashion.
%   
%   See also tsig/subsref, tsig/vertcat, tsig/horzcat.

if nargin<2, K=inf; end % no upper bound

Nc = nchan(T);
Nk = length(K);
if Nc==1,
    error(numericTest(K,'scalar/posint/nonnan', 'Second argument is '));
else,
    error(numericTest(K,'posint/nonnan', 'Second argument is '));
    if ~isequal(1,Nk) && ~isequal(Nc,Nk),
        error('Incompatible number of channels between tsig object and second argument.');
    end
end

% find indices using find, and derive time values from them
dT = dt(T); % time spacing in ms
t0 = onset(T); % onsets in ms per channel
K = SameSize(K,1:Nc);
for ichan=1:Nk,
    ihit = find(T.Waveform{ichan},K(ichan),varargin{:}); % indices
    t{ichan} = t0(ichan)+dT*ihit(:);
end


