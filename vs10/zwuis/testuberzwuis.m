function [isokay, dpfail] = testuberzwuis(DPorder, freq);
% testuberzwuis - test whether freqs array meets uberzwuis criteria
%    testuberzwuis(Freq, DPorder) tests whether the DPs of order DPorder 
%    evoked by frequencies Freq are all different. A logical value is
%    returned indicating the success of the test. DPorder may be an array, 
%    in which case the test is performed not only within the respective 
%    orders (elements of DPorder), but also across orders.
%   
%    [isokay, DPfail] = testuberzwuis(Freq, DPorder) also returns a matrix
%    DPfail whise rows contain the weights of non-unique DPs. E.g., if
%    DPfail(k,:) equals [0 2 0 1 0 0], this means that 2*Freq(2)+1*Freq(4)
%    is a non-unique DP. More generally, the array DPfail*Freq(:) contains
%    those DP freqs that are non-unique.
%
%    See also uberzwuis, oldDPmatrix.

N = numel(freq); % # comps
M = oldDPmatrix(N,DPorder);
DPfreq = M*freq(:);

[DPsort, isort] = sort(DPfreq);
isame = find(diff(DPsort)==0);
isamesort = unique([isame; isame+1]); % indices of non-unique elements in DPsort 
isame = isort(isamesort); % indices of non-unique elements in (unsorted) DPfreq ...
dpfail = M(isame,:); % ... and thus of the corresponding rows in M
isokay = isempty(isame);



