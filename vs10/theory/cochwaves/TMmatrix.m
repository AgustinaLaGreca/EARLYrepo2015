function TM = TMmatrix(N,M,Di)
% TMmatrix - potential energy matrix of transmission line
%   TMmatrix(N,M,Di) returns NxN matrix describing an elastic interaction
%   between neighboring points on a closed transmission line. M is the 
%   order of the interaction, corresponding to the order of DIFF:
%      M=0 yields identity matrix ("elastic base")
%      M=1 sums the differences between adjacent points ("tension")
%      M=2 sums the diff of the diff of adjacent points ("bending stiffness")
%   Di is the distance between the coupled neighbors. Default Di = 1.
%
%   See also DIFF.

Di = arginDefaults('Di',1);
% there will no doubt be an elegant way of doing this, but let's be
% practical for one time.

switch M
    case 0,
        TM = eye(N);
    case 1,
        TM = zeros(N);
        mm = [1 -1; -1 1];
        for ii=1:N,
            iin = 1+rem(ii+Di-1,N); % neighbor on the circle
            TM([ii iin], [ii iin]) = TM([ii iin], [ii iin]) + mm;
        end
    case 2,
        TM = zeros(N);
        mm = [1 -2 1; -2 4 -2; 1 -2 1];
        for ii=1:N,
            iin = 1+rem(ii+Di-1,N); % neighbor on the circle
            iinn = 1+rem(ii+2*Di-1,N); % next neighbor on the circle
            TM([ii iin iinn], [ii iin iinn]) = TM([ii iin iinn], [ii iin iinn]) + mm;
        end
    otherwise,
        error('M must be 0,1, or 2');
end









