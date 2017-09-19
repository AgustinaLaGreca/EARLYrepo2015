function [best_freq, DPcount, DP, Ndp] = uberzwuis(N, Nit, Fmax, DFmin);
% uberzwuis - array of freqs with unique 3rd order distortions
if nargin<4, DFmin = 1; end
if (N-1)*DFmin > Fmax, 
    error('Too many widely spaced components to fit in Fmax.');
end

M = local_DP3mat(N);
Ndp = size(M,1); % # pos DPs
MinDPcount = -inf;
best_freq = nan;
for ii=1:Nit,
    freq = 1+floor(Fmax*rand(N-1,1)); % N freq except the first one, which vanishes by definition
    if any(diff([0; freq])<DFmin), continue; end % 
    DP = M*freq; % all positive-weighted 3rd order DPs
    DPcount = length(unique(DP));
    if DPcount>MinDPcount,
        best_freq = sort(freq);
        MinDPcount = DPcount;
    end
    if isequal(Ndp, DPcount), % all DPs are unique, quite searching 
        break;
    end
end
DP = (M*best_freq)';
best_freq = [0 best_freq'];
DPcount = MinDPcount;


%----------------------
function M = local_DP3mat(N);
% all choices of weights out of N, totalling 3
M = zeros(N^3,N);
icount = 0;
for i1=1:N,
    for i2=1:N,
        for i3=1:N,
            icount = icount + 1;
            M(icount, i1)  = M(icount, i1)  + 1;
            M(icount, i2)  = M(icount, i2)  + 1;
            M(icount, i3)  = M(icount, i3)  + 1;
        end
    end
end
M = unique(M,'rows');
M = M(:,2:end); % first col weights zero freq so may be omitted






