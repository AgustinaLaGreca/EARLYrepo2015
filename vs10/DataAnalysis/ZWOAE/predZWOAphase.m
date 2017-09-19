function predZWOAphase(F1, F2, delay1, delay2, tauAGC, wholeCycle);
% predZWOAphase - predict phases of ZWOA data
%     syntax:
%     predZWOAphase(F1, F2, delay1, delay2, tauAGC);
%
%     Freq in Hz; delays in us.
%     F1 is zwuis group; F2 is lonely component.
%     Delays are delays toward/from DP generation site

if nargin<6, wholeCycle=0; end
Nzw = length(F1);

% convert freqs to kHz
F1 = F1/1e3;
F2 = F2/1e3;

% linearly extrapolate delay(freq) relation
Pdelay = polyfit([mean(F1) F2], [delay1 delay2],1);

% compute phases at DP site
phase1 = local_delay(0, F1, Pdelay);
phase2 = local_delay(0, F2, Pdelay);

% put all primary freqs & phases in col vectors
FRprim = [F1(:); F2];
PHprim = [phase1(:); phase2];

% aux matrices
Md = -EvaDifMat(Nzw); % differences across Nzw components; if F1 are excending, Md*F1 are positive
Ms = [abs(Md); 2*eye(Nzw)];  % sums across Nzw components, including same-term sums
% DP matrices themselves
Mat_near = [Ms, -ones(size(Ms,1),1)]; % F1(i)+F1(j)-F2 cmp
Mat_far =  [-eye(Nzw), 2*ones(Nzw,1)]; % 2*F2-F1 cmp
Mat_slo = [Md, ones(size(Md,1),1)]; % F2 +F1b-F1a components (upper suppr band)
Mat_shi = [-Md, ones(size(Md,1),1)]; % F2 -F1b+F1a components (lower suppr band)

% compute DP freqs & DP phases, including backward traveling delays
[FRnear, PHnear] = local_distort(Mat_near, FRprim, PHprim, Pdelay, wholeCycle);
[FRfar, PHfar] = local_distort(Mat_far, FRprim, PHprim, Pdelay, wholeCycle);
[FRslo, PHslo] = local_distort(Mat_slo, FRprim, PHprim, Pdelay, wholeCycle);
[FRshi, PHshi] = local_distort(Mat_shi, FRprim, PHprim, Pdelay, wholeCycle);


%======PLOT
%[fig, PS] = WASdagFig;
PS = ZWOAEplotStruct

ha1 = subplot(2,1,1);
fr = linspace(0, max(F1+F2), 1e3);
plot(fr, polyval(Pdelay,fr), 'color', 0.8*[1 1 0.9]);
grid on;
ylabel('Travel delay (\mus)')

ha2 = subplot(2,1,2);
plot(FRnear, PHnear, PS.DPnear)
xplot(FRfar, PHfar, PS.DPfar)
xplot(FRslo, PHslo, PS.DPsup)
xplot(FRshi, PHshi, PS.DPsup)
grid on
local_GDwrite(FRnear, PHnear);
local_GDwrite(FRfar, PHfar);
local_GDwrite(FRslo, PHslo);
local_GDwrite(FRshi, PHshi);
xlabel('Frequency (kHz)')
ylabel('Phase (Cycle)');

set(gcf,'position', [667 33 560 689])
XL = xlim;
axes(ha1); xlim(XL);
set(gcf, 'InvertHardcopy', 'off');

%===========
function phase = local_delay(phase, freq, Pdelay, wholeCycle);
if nargin<4, , wholeCycle= 0.5; end
delay = polyval(Pdelay, freq);
phase = phase - 1e-3*delay.*freq; % factor 1e-3: us->ms
% unwrapping is nontrivial - freqs are not necessarily sorted, ...
% ... and should be kept as is.
[dum isort] = sort(freq);
phase = cunwrap(phase(isort)); % temporarilyt sort the phases according to freq
phase(isort) = phase; % restore order
%Mph = mean(phase), wholeCycle
phase = phase+wholeCycle-round(mean(phase)-wholeCycle); % project to cycle around wholeCycle

function  [FR_dp, PH_dp] = local_distort(Matrix, FRprim, PHprim, Pdelay, wholeCycle);
% compute freq & phases resulting from delays.
FR_dp = Matrix*FRprim; 
PH_dp = Matrix*PHprim;
PH_dp = local_delay(PH_dp, FR_dp, Pdelay, wholeCycle);

function local_GDwrite(FRnear, PHnear);
p = polyfit(FRnear, PHnear,1);
GD = round(-1e3*p(1));
ht = text(mean(FRnear), 0.1+mean(PHnear), [num2str(GD) ' \mus']);




