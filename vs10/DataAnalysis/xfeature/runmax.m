function [Y, Imax] = runmax(X,Nwin);
% runmax - running maximum
%   Y = runmax(X, N) returns the maximima of a N-sample-long sliding 
%   window applied to array X:
%
%       Y(K) = max(X(K-M:K+M)) where N = 2*M+1.
%
%   [Y, iPeak] = runmax(X, N) also returns the indices of samples iPeak
%   at X and Y coincide, i.e., iPeak = find(X==Y).
%
%   See also Smoothen.

isrow = size(X,2)>1;
X = X(:);
% The brute force method using convmtx is fast, but easily runs into memory
% problems. Solve this by chopping X into smaller chunks
Nmax = 1e5; % max # samples in convmtx
%Nmax = 5;


Nsam = numel(X);
M = round((Nwin-1)/2);
Nwin = 2*M+1;

NsamMtx = Nwin*(Nsam+Nwin-1); % # elements in convmtx
Nchunk = ceil(NsamMtx/Nmax); % number of chunks needed to stay below Nmax
NsamChunk = floor(Nsam/Nchunk)+Nwin;
% compute start & end indices of chunks
i0 = 1:NsamChunk:Nsam;
i0 = i0(i0<=Nsam);
i1 = i0+NsamChunk+Nwin-2;
i1 = min(i1,Nsam);
Nchunk = numel(i0);
%Nsam, NsamChunk, i0, i1

Y = []; Imax = [];
Xmin = min(X(:));
X = X-Xmin+1; % ensure positiveness of x, otherwise the zero entries of convmtx below can creep in
for ichunk=1:Nchunk,
    x = X(i0(ichunk):i1(ichunk));
    %mm = convmtx(x,Nwin)
    y = max(convmtx(x,Nwin),[], 2); 
    y = y(1+M:end-M);
    if ichunk>1, % chunk partially overlaps w previous chunk
        Noverlap = min(Nwin-1, numel(y));
        Y(end-Noverlap+1:end) = max(Y(end-Noverlap+1:end), y(1:Noverlap)); % in the overlap region take max of the two
        y = y(Nwin:end); 
    end
    Y = [Y; y];
end
if nargout>1,
    Imax = find(X==Y);
end
Y = Y + Xmin-1; 
if isrow,
    Y = Y.'; 
    Imax = Imax.';
end


