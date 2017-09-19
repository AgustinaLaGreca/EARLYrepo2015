function M = oldDPmatrix(N, DPorder);
% oldDPmatrix - matrix summarizing all DPs of a given order
%    M=DPmatrix(Ncomp, DPorder) returns a matrix M whose rows contain
%    positive integer weights of Ncomp freq components adding up to DPorder.
%    If Freq is a column vector containing frequencies, then M*Freq
%    contains all DPs of order DPorder having positive weights.
%    This is used to check unicity of DPs of a given order, including DPs
%    with negative weights. DPorder may not exceed 5. If DPorder is an
%    array, M combines DPs of different orders.
%
%    See also uberzwuis, testuberzwuis.

if any(DPorder>5) || any(DPorder<1),
    error('DPorder may not exceed 5 and must be positive.')
end

DPorder = unique(DPorder);
M = [];
for ii=1:numel(DPorder),
    fcn = ['local_DP' num2str(DPorder(ii)) 'mat'];
    M = [M; feval(fcn,N)];
end

%====locals doing the real work=============================
function M = local_DP1mat(N);
M = eye(N);

function M = local_DP2mat(N);
% all choices of weights out of N, totalling 3
M = zeros(N^2,N);
icount = 0;
for i1=1:N,
    for i2=1:N,
        icount = icount + 1;
        M(icount, i1)  = M(icount, i1)  + 1;
        M(icount, i2)  = M(icount, i2)  + 1;
    end
end
M = unique(M,'rows');

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

function M = local_DP4mat(N);
% all choices of weights out of N, totalling 4
M = zeros(N^4,N);
icount = 0;
for i1=1:N,
    for i2=1:N,
        for i3=1:N,
            for i4=1:N,
                icount = icount + 1;
                M(icount, i1)  = M(icount, i1)  + 1;
                M(icount, i2)  = M(icount, i2)  + 1;
                M(icount, i3)  = M(icount, i3)  + 1;
                M(icount, i4)  = M(icount, i4)  + 1;
            end
        end
    end
end
M = unique(M,'rows');

function M = local_DP5mat(N);
% all choices of weights out of N, totalling 5
M = zeros(N^5,N);
icount = 0;
for i1=1:N,
    for i2=1:N,
        for i3=1:N,
            for i4=1:N,
                for i5=1:N,
                    icount = icount + 1;
                    M(icount, i1)  = M(icount, i1)  + 1;
                    M(icount, i2)  = M(icount, i2)  + 1;
                    M(icount, i3)  = M(icount, i3)  + 1;
                    M(icount, i4)  = M(icount, i4)  + 1;
                    M(icount, i5)  = M(icount, i5)  + 1;
                end
            end
        end
    end
end
M = unique(M,'rows');


