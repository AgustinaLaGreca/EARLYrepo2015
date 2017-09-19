function D = DifMat(N);
% EvaDifMat - matrix for evaluating differences across vector elements

D = zeros(N*(N-1)/2,N);
id = 0; % row counters for D 
for ii=1:N,
    for jj=ii+1:N,
        id = id+1;
        D(id, [ii jj]) = [1 -1];
    end
end


