function [M,combifactor, dbdown] = DPmatrix(Nprim, Norder, SumWeight);
% DPmatrix - weight matrix for the evaluation of distortion products
%   M = DPmatrix(N, Norder) returns a N-column weight matrix M.
%   Let F be the column matrix representing the frequencies F(1)...F(N)
%   of a tone complex, then the totality of DPs of order Norder is given by
%   Fdp = M*F. The rows of M are all different, and pairs of rows are never
%   each other's negations. This excludes the coexistence of weight rows
%   that produce both the positive and the negative frequencies of a given
%   DP component.
%
%   Note that the same weight matrix M can also be used to compute the DP
%   phases from the phases of the primaries.
%
%   If Norder is an array, M combines the different orders by vertical
%   concatenation:
%      M = [DPmatrix(N,Norder(1)); DPmatrix(N,Norder(2)); ...]
%
%   M = DPmatrix(N, Norder, S) restricts M to those rows whose sum equals S.
%
%   [M, CombiFactor, dBdown] = DPmatrix(...) also returns the combinatorial
%   factor CombiFactor and relative levels of the DPs, dBdown (see DPfreqs).
%
%   See also DPfreqs, uberzwuis, ZWOAEmatrices, partitions.

if nargin<3, SumWeight=[]; end

% check for cached values
CFN = mfilename;
CPRM = {Nprim, Norder, SumWeight};
S = getcache(CFN, CPRM);
if ~isempty(S),
    [M,combifactor, dbdown] = deal(S.M, S.combifactor, S.dbdown);
    return;
end

if numel(Norder)>1, % combine orders using recursive call
    [M,combifactor, dbdown] = deal([]);
    for iorder=1:numel(Norder),
        [m,c,d] = DPmatrix(Nprim, Norder(iorder), SumWeight);
        M = [M; m];
        combifactor = [combifactor; c];
        dbdown = [dbdown; d];
    end
    putcache(CFN, 1e3, CPRM, CollectInStruct(M,combifactor, dbdown));
    return;
end

% ======single Norder from here========
M = [];
P = partitions(Norder,Nprim);
Npart = size(P,1);
for ipart=1:Npart,
    p = P(ipart,:);
    p=localSigFlipper(p);
    M = [M; p];
end
% restrict M by applying SumWeight requirement if specified
if ~isempty(SumWeight),
    M = M(sum(M,2)==SumWeight,:);
end

% weed out any "negated pairs" (see help text)
SW = sum(M,2);
M0 = M(SW==0,:); % isolate all zero-sumweight rows
M = M(SW>0,:); % this eliminates the easy ones having nonzero sumweight. 
% Now eliminate the zero-sumweight pairs by brute force.
N0 = size(M0,1);
iweed = []; % death list
for ii=1:N0,
    minRow = -M0(ii,:);
    for jj=ii+1:N0,
        if isequal(minRow, M0(jj,:)),
            iweed = [iweed, jj];
        end
    end
end
M0(iweed,:) = []; 
M = [M0; M]; % recombine with the nonzero=sumweight ones

N = size(M,1); % # DPs
[combifactor, dbdown] = deal(zeros(N,1));
for ii=1:N,
    [combifactor(ii), dbdown(ii)]= local_combi(M(ii,:));
end
putcache(CFN, 1e3, CPRM, CollectInStruct(M,combifactor, dbdown));
%=====================

function p=localSigFlipper(p);
% matrix containing all +- variants of p(1) to p(n) just once.
M = numel(p);
if (M==1) && ~isequal(0,p),
    p = [p; -p];
elseif M>1, % recursion
    q = localSigFlipper(p(2:end));
    NN = size(q,1);
    p1 = ones(NN,1)*p(1);
    if p(1)==0,
        p = [p1 q];
    else,
        p = [p1 q; -p1 q];
    end
end

%=====================
function [f, dBdown] = local_combi(w);
w = abs(w);
w = w(w>0);
dd = prod(factorial(w));
f = factorial(sum(w))/dd;
dBdown = -A2dB(dd);









