function rmin = rayleigh_rmin(N, pval);
% rayleigh_rmin - minimum vector strength meeting pval given # events
% rayleigh_rmin(N, pval) returns the vector strength just meeting the
% criterion p<=pval for a given number N of events, i.e., the r for which 
%    rayleigh_pval(r,N) equals pval.
%  
%    Default N, pval is 0.001. N and N, pval may be arrays of compatible
%    size as described in SameSize.
%
%   See also rayleigh_pval, rayleigh_test.

[pval] = arginDefaults('pval', 0.001);

[N, pval] = SameSize(N, pval);

if numel(N)>1,
    rmin = nan(size(N));
    for ii=1:numel(N),
        rmin(ii) = rayleigh_rmin(N(ii), pval(ii));
    end
    return;
end

%=====single N and pval from here=======
pfcn = @(r)(rayleigh_pval(r,N)-pval);
if prod(pfcn([0 1]))>=0, 
    rmin=nan;
else,
    rmin = fzero(@(r)(rayleigh_pval(r,N)-pval),[0 1]);
end





