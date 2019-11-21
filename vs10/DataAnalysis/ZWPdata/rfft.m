function X = rfft(x)
% rfft - fourier transform of real waveform
%    X = rfft(x) returns the positive-frequency components of the discrete
%    fourier transform of the real-valued array x, normalized such that
%                     N
%       x(n) = real( sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N) )
%                    k=1
%
%    Here, N = length(X) = (M+2)/2,   even M = length(x)
%                        = (M+1)/2,   odd M.
%
%    The logic of this normalization is that A(k)=abs(X(k)) is the linear
%    amplitude of the kth frequency component and that ph0 = angle(X(k)) is
%    its phase: 
%
%        comp(t) = A(k)*cos(omega(k)*t+ph0(k))
%
%    with omega(k) = 2*pi*(k-1)*/Duration.
%    Note that the inverse transform cannot be unambigously formulated
%    because numel(x) cannot be inferred from numel(X).
%
%    If x is a matrix, rfft operates on its columns: X(:,m) = rfft(:,m).
%
%    See also fft, ifft.

if ~isreal(x),
    error('Input array x must be real.');
end
Sz = size(x);
[x, isRow] = TempColumnize(x);
Nsam = size(x,1);
X = fft(x)*(2/Nsam); % factor 2 anticipates the real() operation mentioned in the help text
N = floor((Nsam+2)/2); % positive freqs only.
X(N+1:end,:)=[];
X(1,:) = X(1,:)/2; % undo doubling of DC component(s)
if ~rem(Nsam,2), % force Nyquist cmp to be real
    X(end,:) = real(X(end,:));
end
if isRow,
    X = X.';
end








