function [heq, Ieta2, ST] = BMeq_height(lambda, yBM, dBM, ISL_BM, BM_OSL, L2, N);
% BMeq_height - Steele & Taber's equivalent height of the fluid pushed by the BM
%   BMeq_height(lambda, yBM, dBM, ISL_BM, BM_OSL, L2, N)
%   Specify all input args (except dBM, which will be normalized) in um.
%   Default values indicated below in [brackets].
%     lambda: longitudinal wavelength
%        yBM: radial positions for the BM displacement profile
%        dBM: radial BM displacement profile (shape fcn). Will be normalized.
%     ISL_BM: radial distance from inner wall of scala to inner end of BM [100 um]
%     BM_OSL: radial distance from outer end of BM to outer wall of scala [50 um]
%         L2: height of scala [300 um]
%          N: number of points in spatial spectrum of z=0 profile [2048]
%
%    BMeq_height returns the equivalent height in um of the layer of fluid moved
%    by the BM displacement.
%
%    [h_eq, Ieta2, ST] = BMeq_height(...) also returns the "eta integral" Ieta2, 
%    which is the squared normalized BM displacement integrated over the 
%    radial width (this integral occurs in Steele & Taber Eqs  13, 14b,
%    16b), and a struct ST containing the m_j (radial and transverse
%    wavenumbers, eq. 8) and B_j (coeffients in the development of the
%    velocity potential, eqs. 7 & 9). The rows of ST.m and ST.B correspond
%    to the elements of lambda; the subscript j in S&T's notation runs
%    along the columns.

[ISL_BM, BM_OSL, L2, N] = arginDefaults('ISL_BM/BM_OSL/L2/N', 100, 50, 300, 2048);

sLambda = size(lambda);
lambda = lambda(:);

yBM = yBM(:).' - yBM(1); % row vector; left BM edge ~ y=0 by convention
dBM = dBM(:).'/max(dBM); % row vector; normalize to unity @ peak
mindy = min(diff(yBM));
yBMinner = mindy/1e6; % just inside specified BM pos range
yBMouter = yBM(end)+mindy/1e6; % just outside specified BM pos range
Y = [-ISL_BM, yBMinner, yBM, yBMouter, yBMouter+BM_OSL];
L1 = max(Y)-min(Y);
Eta = [0 0 dBM 0 0]; % radial shape function, including bony rims neighboring the BM

% resample Eta to get regular spacing needed for spatial FFT.
nY = linspace(min(Y), max(Y), N);
dY = diff(nY(1:2));
Eta = interp1(Y,Eta, nY);
% evaluate Ieta2
Ieta2 = dY*sum(Eta.^2); % Steele's eta-integral in e.g. Eq 13.

% Work toward the A_j components in Eq 10a.
Eta = [Eta, 0*Eta]; % add zeros to get half-cycle FFT cmps
Csp = fft(Eta); % i.e., Csp(j+1) = Sum_n Eta(n*dY).*exp(-2*pi*i*j*n/(2*N))
% ... = Sum_y Eta(y).*exp(-pi*i*j*y/(N*dY)) 
% ... = Sum_y Eta(y).*exp(-pi*i*j*y/L1), with y = n*dY = 0..2*L1
rCsp = real(Csp); % force cosine spectrum:
%  rCsp(j+1) = real(Sum_y Eta(y).*exp(-pi*i*j*y/L1))
%            = Sum_y Eta(y).*cos(pi*j*y/L1)
dk = 2*pi/(dY*numel(Csp)); % wavenumber spacing of spectrum
rCsp(N:end) = 0; % Wipe out "neg freqs"
% now use Steele and Taber's (1979) notation
j = 0:N-2;
jpos = 1+(j>0);
A(1,j+1) = rCsp(j+1)*dY; % Eq 10a: A(j+1) = Sum_y Eta(y).*cos(pi*j*y/L1)*dY ...
%                             ... ~= Integral dy Eta(y).*cos(pi*j*y/L1) 

% if 1,
%     % Check whether spatial spectrum yields correct profile
%     NDC2 = 1+((1:numel(rCsp))>1); % doubles all comps except DC; needed to correct discarding of "neg freqs"
%     rEta = 2*real(ifft(NDC2.*rCsp)); % factor of 2 compensates doubling of # spectral components ...
%     % ... (padding zeros). Needed because Matlab's ifft uses a 1/N normalization.
%     subplot(2,1,1)
%     dplot(dY, Eta)
%     xdplot(dY, rEta, 'r');
%     subplot(2,1,2);
%     Psp = a2db(abs(Csp));
%     dplot(dk, Psp);
%     xlim([0 0.5]);
% end


%     -> Real(Sum_k dY exp(2*pi*i*k*j/N) eta(k*j/N)) = fft(eta)*dY. 
k = 2*pi./lambda; % Use k for wavenumber, not Steele's unusual notation lambda
[k, A, j, jpos] = SameSize(k, A, j, jpos); % A, j, and jpos are row arrays; k is col array
m = sqrt((j*pi/L1).^2 + k.^2); % Eq 8. 
heq = (1/Ieta2)*sum(jpos.*A.^2./(m*L1.*tanh(m*L2)),2); % Eq 14b. Sum is over dim 2, i.e. over j=0:N-1

heq = reshape(heq, sLambda);
lambda = reshape(lambda, sLambda);
ST = CollectInStruct('-param', lambda, yBM, dBM, ISL_BM, BM_OSL, L1, L2, N, ...
    '-coeff', k, m, A, heq, Ieta2);















