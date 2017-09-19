function W = ZWOAEfit(Gerbil, idataset, DPtype, PHref, plotflag);
% ZWOAEfit -fit 1st-order polynomial to ZWOAE fase data
%
% W = ZWOAEfit(Gerbil, idataset) takes the specified dataset & fits a plane to
% the phase data of the 'supall' emissions group (re. prim. phase from spectrum')
% In this fit, the fixed-stimulus frequency (non-zwuis) is not considered.
% The plane is specified as: PHdp = PHo + beta(1)*(Fa+Fb) + beta(2)*(Fa-Fb)
%
% Multiple idataset are allowed. This will result in output W to be an array of struct.
%
% W = ZWOAEfit(Gerbil, idataset, DPtype) allows specification of the DPOAE
% group to consider: default = 'supall', other options are: 'near',
%                                       'suplo', 'suphi'
% NOTE: for 'far' (& other options in ZWOAEsubspec) the plane cannot be
% fitted. Therefore, those options are not allowed.
%
% W = ZWOAEfit(Gerbil, idataset, DPtype, PHref) allows specification of the
% DPOAE phase reference. Default = 'response'. Other options are 'stimulus'
% or a numeric value (see ZWOAEsubspec.m)
%
% W = ZWOAEfit(Gerbil, idataset, DPtype, PHref, 1) also generates 2-D plots
% of the phase data + model result
%
% W is (an array of) struct, holding the following fields:
%     'expID',  'recID', 'ZwuisLoc', 'DPtype', 'PHref' (the usual...)
%    'BaseFnc'  : cell array of strings, specifying the fitted param. combinations  
%     'Fr'      : frequencies of the DPs [kHz]
%     'Ph'      : phase of the DPs [cycle]
%    'Ph_recalc': phase of the DPs, calculated from the fitted polynomial [cycle]
%     'beta'    : fit coefficients, given in same order as BaseFnc
%     'VA'      : variance accounted for by the fit [%]
%
%
% See also ZWOAEsubspec, ZWOAEfase
%

%--- check on inputs / assign defaults ---
if nargin < 5 | isempty(plotflag), plotflag = 0; end
if nargin < 4 | isempty(PHref), PHref = 'response'; end
if nargin < 3 | isempty(DPtype), DPtype = 'supall'; end
%validity-check of PHref is done by ZWOAEsubspec

%check validity of DPtype
[DPtype, Mess] = keywordMatch(DPtype, {'near' 'supall' 'suplo' 'suphi'});
error(Mess);

%--- use recursion for multiple idataset
if numel(idataset) > 1,
    for ii = 1:numel(idataset),
        W(ii) = ZWOAEfit(Gerbil, idataset(ii), DPtype, PHref, plotflag);
    end
    return
end

%--- single recording from here ---

%get the ZWOAE-data using ZWOAEsubspec
D = ZWOAEsubspec(Gerbil, idataset, DPtype, PHref);

PH_dp = D.PH(D.isortFreq); %DPOAE phases, "unsorted"
W_dp = D.W(D.isortFreq); %weigth factors, "unsorted"
Fprim = [D.Fzwuis(:); D.Fsingle]; % primary freqs in kHz
Nzw = numel(D.Fzwuis); %# of zwuis components

[M, MP] = ZWOAEmatrices(Nzw, ['M' DPtype]);
% MP = index matrix addressing primaries 
 
% Since only 'near' and 'sup...' are allowed, either column 3 (near) or 2 (sup...) is
% constant.
Fa = Fprim(MP(:,1));
ik = 2; %use 2nd column in MP in case of 'near'
if strncmp('sup', DPtype, 3), ik = 3; end; %in case of 'sup...', use 3rd column
Fb = Fprim(MP(:,ik));

BaseFnc = {'1' '(Fa+Fb)' '(Fa-Fb)'}; %these are the functions used in the fit

VDM = GVD(BaseFnc, {'Fa' 'Fb'}, Fa, Fb);
% construct generalized Vandermonde matrix and fit polynomial
% coefficients
[beta, Rdef] = wlinsolve(VDM, PH_dp, W_dp);
if ~isequal(Rdef, 0),
    error(['Problem is rank deficient (' int2str(Rdef) ') & thus not uniquely ''solvable''']);
end
Ph_recalc = VDM*beta;
VA = VarAccount(PH_dp, Ph_recalc, W_dp);

ZwuisLoc = 'below'; if D.Fz_mean > D.Fsingle, ZwuisLoc = 'above'; end

expID = Gerbil;
recID = idataset;
PHref = D.Ref; %retrieve full name for PHref from struct D
Ph = PH_dp; %DP phases
Fr = M*Fprim;

if nargout < 1, 
    %no output, echo results to command line
    disp('---------')
    for ii=1:length(beta),
        disp([sprintf('%8.3f', beta(ii)) ' * ' BaseFnc{ii}]);
    end
    disp('---------')
    disp(['Var. account: ' num2str(round(10*VA)/10) '%']);
else
    W = CollectInStruct(expID, recID, ZwuisLoc, DPtype, PHref, '-',...
        BaseFnc, '-', Fr, Ph, Ph_recalc, beta, VA);
end


%--- plot ---
%generate 2-D plots for the phase data

if plotflag == 1,
    figure;
    set(gcf,'units', 'normalized', 'position', [0.153 0.325 0.522 0.589]);
    Fafine = linspace(min(Fa), max(Fa), 100);
    Fbfine = linspace(min(Fb), max(Fb), 100);
    [mn mx] = minmax([Fa(:); Fb(:)]);
    mn = floor(10*mn)/10; mx = ceil(10*mx)/10;
    [f1grid, f2grid, PH_2Dinterp] = griddata(Fa, Fb, Ph, Fafine.', Fbfine);
    [f1grid, f2grid, PH_2Drecalc] = griddata(Fa, Fb, Ph_recalc, Fafine.', Fbfine);
    [Phmin, Phmax] = minmax([PH_2Dinterp(:);  PH_2Drecalc(:)]);
    Phmin = floor(Phmin*10)/10; Phmax = ceil(Phmax*10)/10;

    contInt = 0.1; %distance between contours (cycle)
    phLevels = Phmin:contInt:Phmax;

    subplot(2,2,1); cla;
    contourf(f1grid, f2grid, PH_2Dinterp, phLevels);
    xlabel('F_a [kHz]'); ylabel('F_b [kHz]'); 
    title('Data');
    axis square, axis([mn mx mn mx]);
    caxis([Phmin Phmax]);

    subplot(2,2,2); cla;
    contourf(f1grid, f2grid, PH_2Drecalc, phLevels);
    xlabel('F_a [kHz]'); ylabel('F_b [kHz]');
    title('Model');
    axis square, axis([mn mx mn mx]);
    caxis([Phmin Phmax]);

    %fake axes to place text with info
    subplot(2,2,3); cla, box off; 
    xlim([0 1]); ylim([0 1]);
    dfc = get(0,'defaultfigurecolor'); %default figure color to "hide" axes
    set(gca,'Xtick',[],'Ytick',[],'color',dfc, 'xcolor', dfc, 'ycolor', dfc);

    %add text stating: # omitted points, group delay, variance accounted
    % and reference phase
    
    if isnumeric(PHref), PHref = num2str(PHref); end
    Infostr = {['expID: ' int2str(Gerbil) '; recID: ' int2str(idataset) '.']};
    Infostr = [Infostr; {['DPtype: ' DPtype '; PHref: ' PHref '.']}];
    Infostr = [Infostr; {['ZwuisLoc: ' ZwuisLoc '.']}];
    Infostr = [Infostr; {['contour interval: ' num2str(contInt) ' cycle.']}];
    Infostr = [Infostr; {'--------------------------------------------------'}];
    Infostr = [Infostr; {'Fitted model:'}];
    Infostr = [Infostr; {' '}];
    for ii = 1:numel(BaseFnc),
        Infostr = [Infostr; {[sprintf('%8.3f', beta(ii)) ' * ' BaseFnc{ii}]}];
    end
    Infostr = [Infostr; {' '}];
    Infostr = [Infostr; {['Var. Acc.: ' num2str(round(VA*10)/10) ' %.']}];
    
    text(.05,.95, Infostr,'verticalalignment','top');
    
    subplot(2,2,4); cla;
    xlim([0 1]); ylim([0 1]);
    axis square; set(gca,'Xtick',[],'Ytick',[]);
    xplot([0.2 0.8],[0.2 0.8],'k');
    text(0.5, 0.9,'Fa + Fb = {\it constant}');
    xplot([0.8 0.2],[0.2 0.8],'r');
    text(0.5, 0.1,'Fa - Fb = {\it constant}','color','r');
   
    figure(gcf);
end



if 1, return; end
%=========================================================================
%===============    Old function from here     ===========================
%=========================================================================

% ZWOAEfit - fit ZWOAE data 
%   syntax:
%    ZWOAEfit(Gerbil, idataset, DPtype, N1, N2);
%    DPtype is one of Near, Far, Suplo Suphi, Supall, All.

if nargin<3, DPtype = 'near'; end
if nargin<4, N1=2; end % default order of tau fit
if nargin<4, N2=2; end % default order of poly fit

[DPtype, Mess] = keywordMatch(DPtype, {'near' 'far' 'suplo' 'suphi' 'supall' 'all'});
error(Mess);
PS = ZWOAEplotStruct; % plotting defaults
PlotSpec = PS.(['DP' DPtype]);
DPtype = ['M' DPtype]; % naming according to ZWOAEmatrices

% get stimpars & spectral data
D = getZWOAEdata(Gerbil, idataset);
F1 = D.stimpars.F1/1e3; % kHz
F2 = D.stimpars.F2/1e3; % kHz
if length(F2)>1, [F1, F2] = swap(F1, F2); end % per convention zwuis group =F1
Fprim = [F1(:); F2]; % primary freqs in kHz
Nzw = length(F1); % # zwuis components in F1 group
[df, MG PH Nsam Pnf] = ZWOAEspec(D); % magn & phase spectra, noise floor

% get correct DP matrix, compute DP freqs & get spectrum @ these freqs
[SM, SF] = ZWOAEmatrices(Nzw);
M = SM.(DPtype); % select requested DP type
Freq_dp = M*Fprim;
MG_dp = local_SpecComp(df, MG, Freq_dp); % magn in dB of DPs
Noise_dp = polyval(Pnf, Freq_dp); % estimated noise floor
% dplot(df, MG,PS.Spec);
% xplot(Freq_dp, MG_dp, PlotSpec);
% xplot(Freq_dp, Noise_dp, 'k-');
SNR_dp = MG_dp - Noise_dp; % S/N ratio of DP data
W_dp = SNR2W(SNR_dp); % weight factors of spectral data @ DP freqs

% get DP phases 
PH_prim = local_SpecComp(df, PH, Fprim); % primary phases
PH_dp = local_SpecComp(df, PH, Freq_dp)-M*PH_prim; % DP phases re prim phases
PH_dp = ucunwrap(PH_dp, Freq_dp); % unwrap along freq axis

% construct "vandermonde" matrix: m-th column is f^m+g^m-h^m + (f+g-h)^m,
% where f,g,h are the 3 primary frequencies behind the DP at freq f+g-h.
iF = SF.(DPtype); % freq indices evoking the DPs of requested DPtype
Freq = Fprim(iF); % the j-th row contains the 3 freqs evoking DP @ Freq_dp(j)
Ndp = size(Freq,1); % # DPs
VDM = zeros(Ndp, N1+1); % N(1) is the order of the polynomial; m=0..N
for m=0:N1,
    VDM(:,m+1) = Freq(:,1).^m + Freq(:,2).^m - Freq(:,3).^m + Freq_dp.^m;
end
VDM = fliplr(VDM); % now highest order is leftmost column, agreeing with ..
% .. Matlab polynomial conventions

% the model is: Phase(f,g,h) = Sum_m alpha(m)*(f^m+g^m-h^m + (f+g-h)^m)
% with unkown coefficients alpha to be optimized.
% Using the Vandermonde matrix, this becomes an overdetermined matrix
% equation: Phase = VDM*alpha, which is solved using wlinsolve.
alpha = wlinsolve(VDM, PH_dp, W_dp);
ph_recalc = VDM*alpha;

% Second type of "model": 1st- or 2nd-order polynomial in the primary freqs that is
% symmetrical in the two frequencies Fa and Fb occurring in the expression 
% Fdp = Fa+Fb-Fc.
MP = SF.(DPtype); % index matrix addressing primaries 
Fa=Fprim(MP(:,1)); Fb=Fprim(MP(:,2)); Fc=Fprim(MP(:,3));
% compute the vandermonde matrix
if N2>2, error('Fitting orders>2 not yet implemented'); end
BaseFnc = {'1' '(Fa+Fb+Fc)/3' '(Fa+Fb-Fc)' '(Fa+Fb-mean(Fa+Fb)).^2' '(Fa-Fb-mean(Fa-Fb)).^2' '(Fc-mean(Fc)).^2' '(Fa+Fb-mean(Fa+Fb)).*(Fc-mean(Fc))'};
%           1  2             3    4                             5                        6                  7   
if N2==1, BaseFnc = BaseFnc(1:3); end % linear terms only

% Note that within a single DPtype only 1 or 2 of the 3 primaries
% frequencies are varied. It is not useful to fit along a fixed parameter.
% So we remove the columns of the VanderMonde that are linearly dependent
% on other columns, thereby fixing any rank deficiency in VDM2.
iremove = []; % coulmns of VDM2 to be removed
aconst = isconstant(Fa);
bconst = isconstant(Fb);
cconst = isconstant(Fc);
Nfvar = 3-aconst-bconst-cconst;
if cconst, iremove = [iremove 2 6 7]; end
if aconst && bconst, iremove = [iremove 2 4 5 7]; end
if aconst || bconst, iremove = [iremove 5]; end
iremove = iremove(iremove<=length(BaseFnc));
BaseFnc(iremove) = [];
VDM2 = GVD(BaseFnc, {'Fa' 'Fb' 'Fc'}, Fa, Fb, Fc);
% construct generalized Vandermonde matrix and fit polynomial
% coefficients
[beta, Rdef] = wlinsolve(VDM2, PH_dp, W_dp); Rdef
disp('---------')
for ii=1:length(beta), disp([sprintf('%8.3f', beta(ii)) ' * ' BaseFnc{ii}]); end
disp('---------')
ph_recalc2 = VDM2*beta;
VA2 = VarAccount(PH_dp, ph_recalc2, W_dp)
%===========PLOT===============================
set(gcf,'units', 'normalized', 'position', [0.43125 0.43167 0.51 0.4975])
[Fmin, Fmax] = minmax([Fprim; Freq_dp]);
freq = linspace(Fmin, Fmax, 1e3);
%---- amplitude
subplot(2,2,1);
dplot(df, MG, PS.Spec);
MG_prim = local_SpecComp(df, MG, Fprim);
MG_dp = local_SpecComp(df, MG, Freq_dp);
xplot(Freq_dp, MG_dp, PlotSpec);
xplot(Fprim, MG_prim, PS.Prim);
xlim([min(Freq_dp) max(Freq_dp) ]);
NF = min(polyval(Pnf, Freq_dp));
ylim([NF inf])
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB SPL)')
title([D.EXPname '/' num2str(D.index) ' ' DPtype])
%----phase data & model fit---
subplot(2,2,2);
plot(Freq_dp, PH_dp, PlotSpec);
xplot(Freq_dp, ph_recalc, 'cs', 'markersize', 3);
xplot(Freq_dp, ph_recalc2, 'ko', 'markersize', 3);
xlim([min(Freq_dp) max(Freq_dp) ]);
xlabel('Frequency (kHz)');
ylabel('DP phase');
%---phase(freq) in the model
subplot(2,2,3);
[Fmin, Fmax] = minmax([Fprim; Freq_dp]);
freq = linspace(Fmin, Fmax, 1e3);
phas = polyval(alpha, freq);
plot(freq, phas);
dfr = diff(freq(1:2));
phas_prim = local_SpecComp([dfr Fmin], phas, Fprim);
phas_dp = local_SpecComp([dfr Fmin], phas, Freq_dp);
xplot(Fprim, phas_prim, PS.Prim);
xplot(Freq_dp, phas_dp, PlotSpec);
grid on
xlabel('Frequency (kHz)');
ylabel('Model phase (cycle)');
%---delay(freq) in the model
subplot(2,2,4);
alpha_deriv = PolyDiff(alpha);
delay = -polyval(alpha_deriv, freq);
plot(freq, delay);
delay_prim = local_SpecComp([dfr Fmin], delay, Fprim);
delay_dp = local_SpecComp([dfr Fmin], delay, Freq_dp);
xplot(Fprim, delay_prim, PS.Prim);
xplot(Freq_dp, delay_dp, PlotSpec);
grid on
xlabel('Frequency (kHz)');
ylabel('Model delay (ms)');

% If Nfvar==2, plot phase data in 2D plot
if Nfvar==2,
    Fpr{1} = Fa; Fpr{2} = Fb; Fpr{3} = Fc; 
    AxLabels = {'Fa' 'Fb' 'Fc'};
    for ii=1:3, Fprfine{ii} = linspace(min(Fpr{ii}), max(Fpr{ii}), 100); end
    if aconst,     i0=1; i1=2; i2=3;
    elseif bconst, i0=2; i1=3; i2=1;
    elseif cconst, i0=3; i1=1; i2=2;
    end
    %dsize(Fpr{i1}, Fpr{i2}, PH_dp)
    [f1grid, f2grid, PH_2Dinterp] = griddata(Fpr{i1}, Fpr{i2}, PH_dp, Fprfine{i1}.', Fprfine{i2});
    [f1grid, f2grid, PH_2Drecalc] = griddata(Fpr{i1}, Fpr{i2}, ph_recalc2, Fprfine{i1}.', Fprfine{i2});
    [Phmin, Phmax] = minmax([PH_2Dinterp(:);  PH_2Drecalc(:)]);
    phLevels = linspace(Phmin, Phmax, 10);
    subplot(2,2,3); cla;
    contourf(f1grid, f2grid, PH_2Drecalc, phLevels);
    xlabel(AxLabels{i1}); ylabel(AxLabels{i2}); 
    subplot(2,2,4); cla;
    contourf(f1grid, f2grid, PH_2Dinterp, phLevels);
    xlabel(AxLabels{i1}); ylabel(AxLabels{i2}); 
end
% F1a, F1b, PH_dp_Mat
% dsize(F1a, F1b, PH_dp_Mat);
% figure


%=========LOCALS=======================
function S = local_SpecComp(df, Spec, Freq);
% return components of spectrum S @ requested frequencies;
if length(df)>1, % 2nd element of df is offset
    offset = df(2);
    df = df(1);
else,
    offset = 0;
end
ii = 1+round((Freq-offset)/df); % indices of components in spectrum
S = Spec(ii);






