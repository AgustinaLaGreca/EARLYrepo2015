function t = time(T, k)
% tsig/time - time values of tsig object
%
%   time(T) or T.time is cell array of time [ms] vectors of the respective
%   channels of tsig T.
%
%   time(T,k) or T.time{k}, where k is a single index, returns the time 
%   values of channel k in a column array.
%
%   See tsig, "methodshelp tsig".

if nargin<2, 
    cellOut = 1; % return cell array (see help)
    k = 1:nchan(T); % default: all channels
else,
    cellOut = 0; % return numeric array (see help)
    if ~isscalar(k), error('Index k must be scalar.'); end
    error(numericTest(k, 'posint/nonnan/noninf', 'Index k is '));
end 

dt = 1/fsam(T); % sample period [ms]
N = nsam(T); % vector; N(k) is channel k
for ii=k,
    t{ii} = T.t0(ii)+ dt*(0:N(ii)-1).';
end

if ~cellOut, t = t{k}; end;



