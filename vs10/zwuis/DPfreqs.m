function [S, DPfreq] = DPfreqs(F, Nmax, WeightSum)
% DPfreqs - low-order distortion products from a set of primary frequencies
%   S = DPfreqs(Freq, Nmax) evaluates distortion products arising by the
%   nonlinear interaction of pure tones at frequencies Freq(k). Nmax is the
%   maximum order of the DPs considered. The return argument S is a struct
%   with fields
%          freq: frequency of the DP in the same units as the input Freq
%         order: order of the SP, which equals sum(abs(weight)).
%        weight: integer weighting factors such that freq = dot(weight,freq)
%     sumweight: sum of weights
%   combifactor: combinatorial factor indicating the algebraic multiplicity 
%                of term coresponding to this DP frequency.
%        dbdown: the amplitude factor in dB corresponding to combifactor.
%          mult: multiplicity, i.e. how many times does freq occur in the
%                complete collection of DPs upto order Nmax.
%   Only non-negative DP frequencies are considered.
%
%   DPfreqs(Freq, Nmax, S) restricts the DPs to those whose summed weights
%   equal S. The default is S=[] meaning that no restriction is imposed.
%
%   Note that weights can also be used to compute the DP phase from the
%   phases of the primaries.
%
%   [S, DPfreq]  = DPfreqs(Freq, Nmax) also returns all the DP frequencies
%   in array DPfreq. That is, DPfreq = [S.freq].
%
%   See also DPmatrix, partitions.
WeightSum = arginDefaults('WeightSum', []); % default: no weight restriction
F = F(:);
Nprim = numel(F);
[M, cfac, dbd] = DPmatrix(Nprim,1:Nmax, WeightSum);
DPfreq = M*F; 
Ndp = numel(DPfreq); % # DPs
% correct any negative elements of DPfreq by negating their weights in M
ineg = find(DPfreq<0);
M(ineg,:) = -M(ineg,:);
DPfreq = abs(DPfreq);
% prepare dealing the rows of M over elements of S
Mc = num2cell(M,2); % each cell is a row
[S(1:Ndp,1).freq] = DealElements(DPfreq);
[S(1:Ndp).order] = DealElements(sum(abs(M),2));
[S(1:Ndp).weight] = deal(Mc{:});
[S(1:Ndp).sumweight] = DealElements(sum(M,2));
[S(1:Ndp).combifactor] = DealElements(cfac);
[S(1:Ndp).dbdown] = DealElements(dbd);
% count multiplicities
[S(1:Ndp).mult] = DealElements(ones(1,Ndp)); % default: assume uniqueness
DPrf = deciRound(DPfreq,7); % when comparing freqs, allow for small rounding errors
[dum, Ifirst, J] = unique(DPrf,'first');
[dum, Ilast, J] = unique(DPrf,'last');
imult = Ifirst(Ifirst~=Ilast); % indices of non-unique DPrf elements
%numel(imult)
for ii = imult(:).',
    qmore = DPrf(ii)==DPrf;
    [S(qmore).mult] = deal(sum(qmore));
end
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









