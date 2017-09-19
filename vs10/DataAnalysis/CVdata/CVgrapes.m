% --------------CVgrapes------
if 1, error('CVgrapes is not a script to be run from commandline.'); end

DS = getDS('RG11337'); % ====================================
f1; CVhi = CVprofile(DS(26:33),2,0,[],-55:5:15); 
f2; CVlo = CVprofile(DS(44:51),2,0,[],-55:5:15); 
f3; CVon = CVprofile(DS(34:41),2,0,[],-55:5:15); 
f4; CVspl = CVprofile(DS(17:25),2,0,[],-55:5:15);

% 1st attempt to =normalize. Use 80-dB zwuis
HI = ds2trf(DS(25),2,0); % 80-dB flat zwuis
Bdr = -30:5:20;
f1; CVhi = CVprofile(DS(26:33),2,0,HI,Bdr); 
f2; CVlo = CVprofile(DS(44:51),2,0,HI,Bdr); 
f3; CVon = CVprofile(DS(34:41),2,0,HI,Bdr); 
f4; CVspl = CVprofile(DS(17:25),2,0,HI,Bdr);



DS = getds('RG11342'); %================================
CM = ds2trf(DS(60),1,-timelag(anachan(DS(60),1))); % use CM as reference
apple(DS(15:22),2,0,CM)

f1; CVspl = CVprofile(DS(17:25),2,0,HI,Bdr);

DS = getds('RG11347'); %=================================
CM = ds2trf(DS(34:36),1,0) % note strange jump @ 25 kHz. CM measurements are strangle SPL-insensitive ..
figure; apple(DS(55:63),2,0,CM); % flat profile, SPL varied
structdisp(stimlist(DS(65:73))) % overview of stim cond
figure; % low-side suppr
apple(DS([65 66 68 70 72]),2,1,CM); % SPL of low-side suppr varied Fsup = 6.1 kHz
apple(DS([65 66 68 70 72]),2,2,CM); % SPL of low-side suppr varied Fsup = 6.5 kHz
apple(DS([65 66 68 70 72]),2,3,CM); % SPL of low-side suppr varied Fsup = 7.1 kHz
apple(DS([65 66 68 70 72]),2,4,CM); % SPL of low-side suppr varied Fsup = 7.5 kHz
apple(DS([65 66 68 70 72]),2,5,CM); % SPL of low-side suppr varied Fsup = 8.1 kHz
figure; % high-side suppr
apple(DS([64 67 69 71 73]),2,1,CM); % SPL of high-side suppr varied Fsup = 15.3 kHz
apple(DS([64 67 69 71 73]),2,2,CM); % SPL of high-side suppr varied Fsup = 15.8 kHz
apple(DS([64 67 69 71 73]),2,3,CM); % SPL of high-side suppr varied Fsup = 16.3 kHz
apple(DS([64 67 69 71 73]),2,4,CM); % SPL of high-side suppr varied Fsup = 16.9 kHz
apple(DS([64 67 69 71 73]),2,5,CM); % SPL of high-side suppr varied Fsup = 17.3 kHz
apple(DS([64 67 69 71 73]),2,6,CM); % SPL of high-side suppr varied Fsup = 17.7 kHz
% flat profi


DS = getds('RG11349'); %=================================
ME = ds2trf(DS([16:18]),2,0) % 80/80/65 dB SPL
AC = ds2trf(DS([16:18]),1,0) % 80/80/65 dB SPL
% pure tones; normalized amp plot
magn_phase_plot(DS(96),2, 'DoNormalize', 1);
magn_phase_plot(DS(96),2, 'DoNormalize', 1, 'GainShift', 80);
apple(DS(49:57),2); % wideband version. Small, but systematic gain effect in LF tail?
ilow = [61 62 66 69 73 74 77:79];
imid = [59 63 67 68 72 75];
ihig = [60 64 65 70 71 76];
AA = apple(DS([ilow imid ihig]),2); % all 1-comp supp data
Fprofile = unique([AA.Fprofile])
baseSPL = unique([AA.baseSPL])
SPLjump = unique([AA.SPLjump])
%-
apple(DS(ilow),2,1);  % 6.1-kHz supp
apple(DS(ilow),2,2);  % 6.5-kHz supp
apple(DS(ilow),2,3);  % 7.1-kHz supp
apple(DS(ilow),2,4);  % 7.5-kHz supp
apple(DS(ilow),2,5);  % 8.1-kHz supp
%-
apple(DS(imid),2,1);  % 15.3-kHz supp
apple(DS(imid),2,2);  % 15.8-kHz supp
apple(DS(imid),2,3);  % 16.3-kHz supp
apple(DS(imid),2,4);  % 16.9-kHz supp
apple(DS(imid),2,5);  % 17.3-kHz supp
apple(DS(imid),2,6);  % 17.7-kHz supp
%-
apple(DS(ihig),2,1);  % 18.3-kHz supp
apple(DS(ihig),2,2);  % 18.8-kHz supp
apple(DS(ihig),2,3);  % 19.5-kHz supp
apple(DS(ihig),2,4);  % 19.9-kHz supp
apple(DS(ihig),2,5);  % 20.5-kHz supp
apple(DS(ihig),2,6);  % 20.8-kHz supp
%-
AAref = AA([AA.SPLjump]<24 & ~betwixt([AA.Fprofile], 8e3, 25e3))
refGain = mean(cat(1,AAref.Gain),1)
AAref = apple(DS(50),2);
refGain = AAref.Gain;
% plot(cat(1,AAref.Fprim).', cat(1,AAref.Gain).')
% legend({AAref.LegStr})
% xplot(AAref(1).Fprim, refGain, 'k', 'linewidth',2)

for ii=1:numel(AA),
    AA(ii).Suppression = refGain - AA(ii).Gain;
end

ifreq = 1:17;
fh(1)=figure; fh(2)=figure;
Bdr = -20:10:50;
for ii=1:numel(ifreq),
    figure(fh(1+floor(ii/10)));
    isub = 1+mod(ii-1,9);
    fpro = Fprofile(ifreq(ii));
    AAf = AA([AA.Fprofile] == fpro); 
    AAf = sortaccord(AAf,[AAf.SPLjump]); 
    GG = cat(1,AAf.Gain); 
    SS = cat(1,AAf.Suppression);
    SS(1,1)=-20.1; SS(2,1)=50.1;
    subplot(3,3,isub); plot([AAf.SPLjump], SS);
    [CS, H] = contourf([AAf(1).Fprim]/1e3, [AAf.SPLjump], SS, Bdr);
    xlog125([5 25]);
    set(gca,'xtick',[5:10 12 15 17 20 25])
    grid on
    clabel(CS,H)
    %colorbar;
    title(sprintf('Fsup= %0.1f kHz', fpro/1e3));
end
% check non-coincidence of 2nd-order (sum) distortions w any primary freqs
Fzw = DS(50).Stim.Fzwuis; Fzw = round(100*Fzw); 
SM = diffmatrix(Fzw, -Fzw); SM = unique(SM(:));
any(ismember(SM, Fzw))


DS = getds('RG11352'); %=================================
ME = ds2trf(DS([5 6 7]),2,0) % 80/80/65 dB SPL
ME2 = ds2trf(DS([19 18 7]),2,0) % 50/50/65 dB SPL
CM = ds2trf(DS([18 9 19]),1,0); % 50 dB SPL
% pure tones; normalized amp plot
magn_phase_plot(DS(28),2, 'DoNormalize', 1);
% magn_phase_plot(DS(28),2, 'DoNormalize', 1, 'GainShift', 80);
% ==flat zwuis====
apple(DS(29:37),2); % uncalibrated
apple(DS(29:37),2, 0, ME); % re middle ear @ 80 dB SPL
apple(DS(29:37),2, 0, ME2); % re middle ear @ 50 dB SPL
apple(DS(29:37),2, 0, CM); % re CM @ 50 dB SPL
% note sharp bend in phase plot @ 9.6 kHz when using either ME or CM; to
% view it, use 1.2-ms delay. Residual delays are then 80 us below 9.6 kHz
% and 230 us between 9.6 kHz and the plateau at ~13 kHz.

% ========suppression data========
ilo = [48 53 54 60 65 66 73];
iml = [49 50 56 59 63 71]; % 67, but this one is corrupt
imh = [47 52 55 61 62 68 72];
ihi = [46 51 57 58 64 69 70];
AA = apple(DS([ilo iml imh ihi]),2); % all 1-comp supp data
Fprofile = unique([AA.Fprofile])
baseSPL = unique([AA.baseSPL])
SPLjump = unique([AA.SPLjump])
%-
apple(DS(ilo),2,1,qq);  % 6.1-kHz supp
apple(DS(ilo),2,2,qq);  % 6.5-kHz supp
apple(DS(ilo),2,3,qq);  % 7.1-kHz supp
apple(DS(ilo),2,4,qq);  % 7.5-kHz supp
apple(DS(ilo),2,5,qq);  % 8.1-kHz supp
apple(DS(ilo),2,6,qq);  % 8.1-kHz supp
%-
apple(DS(iml),2,1);  % 15.3-kHz supp
apple(DS(iml),2,2);  % 15.8-kHz supp
apple(DS(iml),2,3);  % 16.3-kHz supp
apple(DS(iml),2,4);  % 16.9-kHz supp
apple(DS(iml),2,5);  % 17.3-kHz supp
apple(DS(iml),2,6);  % 17.7-kHz supp
%-
apple(DS(imh),2,1);  % 15.3-kHz supp
apple(DS(imh),2,2);  % 15.8-kHz supp
apple(DS(imh),2,3);  % 16.3-kHz supp
apple(DS(imh),2,4);  % 16.9-kHz supp
apple(DS(imh),2,5);  % 17.3-kHz supp
apple(DS(imh),2,6);  % 17.7-kHz supp
%-
apple(DS(ihi),2,1);  % 18.3-kHz supp
apple(DS(ihi),2,2);  % 18.8-kHz supp
apple(DS(ihi),2,3);  % 19.5-kHz supp
apple(DS(ihi),2,4);  % 19.9-kHz supp
apple(DS(ihi),2,5);  % 20.5-kHz supp
apple(DS(ihi),2,6);  % 20.8-kHz supp
%-
AAref = AA([AA.SPLjump] == 0);
refGain = mean(cat(1,AAref.Gain),1)

for ii=1:numel(AA),
    AA(ii).Suppression = refGain - AA(ii).Gain;
end

ifreq = 1:24;
fh(1)=figure; fh(2)=figure; fh(3)=figure;
Bdr = -10:10:40;
for ii=1:numel(ifreq),
    figure(fh(1+floor(ii/10)));
    isub = 1+mod(ii-1,9);
    fpro = Fprofile(ifreq(ii));
    AAf = AA([AA.Fprofile] == fpro); 
    AAf = sortaccord(AAf,[AAf.SPLjump]); 
    GG = cat(1,AAf.Gain); 
    SS = cat(1,AAf.Suppression);
    SS(1,1)=-20.1; SS(2,1)=50.1;
    subplot(3,3,isub); plot([AAf.SPLjump], SS);
    [CS, H] = contourf([AAf(1).Fprim]/1e3, [AAf.SPLjump], SS, Bdr);
    xlog125([5 25]);
    set(gca,'xtick',[5:10 12 15 17 20 25])
    grid on
    clabel(CS,H)
    %colorbar;
    title(sprintf('Fsup= %0.1f kHz', fpro/1e3));
end
% check non-coincidence of 2nd-order (sum) distortions w any primary freqs
Fzw = DS(50).Stim.Fzwuis; Fzw = round(100*Fzw); 
SM = diffmatrix(Fzw, -Fzw); SM = unique(SM(:));
any(ismember(SM, Fzw))


DS = getds('RG11353'); %=================================
ME = ds2trf(DS([1 2 3]),2,0) % 80/80/65 dB SPL
ME2 = ds2trf(DS([5 6 7]),2,0) % 50 dB SPL
CM = ds2trf(DS([5 6 7]),1,0); % 50 dB SPL
% pure tones; normalized amp plot
magn_phase_plot(DS(28),2, 'DoNormalize', 1);
% magn_phase_plot(DS(28),2, 'DoNormalize', 1, 'GainShift', 80);
% ==flat zwuis====
apple(DS(29:37),2); % uncalibrated
apple(DS(29:37),2, 0, ME); % re middle ear @ 80 dB SPL
apple(DS(29:37),2, 0, ME2); % re middle ear @ 50 dB SPL
apple(DS(29:37),2, 0, CM); % re CM @ 50 dB SPL
% note sharp bend in phase plot @ 9.6 kHz when using either ME or CM; to
% view it, use 1.2-ms delay. Residual delays are then 80 us below 9.6 kHz
% and 230 us between 9.6 kHz and the plateau at ~13 kHz.

% ========suppression data========
ilo = [38 45 46 51 56 62 63 68 69 72 73];
iml = [43 47 48 53 58 59 66 67 70 71 74]; % 39 also in this series, but has different random seed for stim freqs
imh = [40 44 49 54 55 60 65];
ihi = [41 42 50 52 57 61 64];
AAlo = apple(DS(ilo),2);
AAml = apple(DS(iml),2);
AAmh = apple(DS(imh),2);
AAhi = apple(DS(ihi),2);
AA = [AAlo AAml AAmh AAhi]; % all 1-comp supp data
Fprofile = unique([AA.Fprofile])
baseSPL = unique([AA.baseSPL])
SPLjump = unique([AA.SPLjump])
%-
apple(DS(ilo),2,1);  % 6.1-kHz supp
apple(DS(ilo),2,2);  % 6.5-kHz supp
apple(DS(ilo),2,3);  % 7.1-kHz supp
apple(DS(ilo),2,4);  % 7.5-kHz supp
apple(DS(ilo),2,5);  % 8.1-kHz supp
apple(DS(ilo),2,6);  % 8.1-kHz supp
%-
apple(DS(iml),2,1);  % 15.3-kHz supp
apple(DS(iml),2,2);  % 15.8-kHz supp
apple(DS(iml),2,3);  % 16.3-kHz supp
apple(DS(iml),2,4);  % 16.9-kHz supp
apple(DS(iml),2,5);  % 17.3-kHz supp
apple(DS(iml),2,6);  % 17.7-kHz supp
%-
apple(DS(imh),2,1);  % 15.3-kHz supp
apple(DS(imh),2,2);  % 15.8-kHz supp
apple(DS(imh),2,3);  % 16.3-kHz supp
apple(DS(imh),2,4);  % 16.9-kHz supp
apple(DS(imh),2,5);  % 17.3-kHz supp
apple(DS(imh),2,6);  % 17.7-kHz supp
%-
apple(DS(ihi),2,1);  % 18.3-kHz supp
apple(DS(ihi),2,2);  % 18.8-kHz supp
apple(DS(ihi),2,3);  % 19.5-kHz supp
apple(DS(ihi),2,4);  % 19.9-kHz supp
apple(DS(ihi),2,5);  % 20.5-kHz supp
apple(DS(ihi),2,6);  % 20.8-kHz supp
%-
AAref = AA([AA.SPLjump] == 0);
refGain = mean(cat(1,AAref.Gain),1)

for ii=1:numel(AA),
    AA(ii).Suppression = refGain - AA(ii).Gain;
end
aa
ifreq = 1:24;
fh(1)=figure; fh(2)=figure; fh(3)=figure;
Bdr = -20:10:40;
for ii=1:numel(ifreq),
    figure(fh(1+floor((ii-1)/9)));
    isub = 1+mod(ii-1,9);
    fpro = Fprofile(ifreq(ii));
    AAf = AA([AA.Fprofile] == fpro); 
    AAf = sortaccord(AAf,[AAf.SPLjump]); 
    GG = cat(1,AAf.Gain); 
    SS = cat(1,AAf.Suppression);
    %SS(1,1)=-20.1; SS(2,1)=50.1;
    subplot(3,3,isub); %plot([AAf.SPLjump], SS);
    [CS, H] = contourf([AAf(1).Fprim]/1e3, [AAf.SPLjump], SS, Bdr);
    %surf([AAf(1).Fprim]/1e3, [AAf.SPLjump], SS, Bdr);
    %view(2)
    caxis(gca, [min(Bdr) max(Bdr)])
    xlog125([5 25]);
    set(gca,'xtick',[5:10 12 15 17 20 25])
    grid on
    clabel(CS,H)
    %colorbar;
    title(sprintf('Fsup= %0.1f kHz', fpro/1e3));
end
% check non-coincidence of 2nd-order (sum) distortions w any primary freqs
Fzw = DS(50).Stim.Fzwuis; Fzw = round(100*Fzw); 
SM = diffmatrix(Fzw, -Fzw); SM = unique(SM(:));
any(ismember(SM, Fzw))

%========================================================
DS = getds('RG11354'); %=================================
ME = ds2trf(DS([1 2 3]),2,0) % 80/80/65 dB SPL
plot(ME);
% pure tones; norm. ampl plot. Not a very sens. prep: compression > 50 dB
magn_phase_plot(DS(29),2, 'DoNormalize', 1); % 0-40 dB
magn_phase_plot(DS(28),2, 'DoNormalize', 1); % 50-80 dB
% flat zwuis
iflat = 30:38;
apple(DS(iflat),2); % uncalibrated
apple(DS(iflat),2, 0, ME); % re middle ear @ 80 dB SPL
% ==suppression conditions, including reference conditions
isupp = 39:77;
Stim = stimparam(DS(isupp));
Lowprofile = unique([Stim.ProfileLowFreq]);
ilo = isupp([Stim.ProfileLowFreq]==Lowprofile(1) & [Stim.ProfileSPLjump]>0);
iml = isupp([Stim.ProfileLowFreq]==Lowprofile(2) & [Stim.ProfileSPLjump]>0);
imh = isupp([Stim.ProfileLowFreq]==Lowprofile(3) & [Stim.ProfileSPLjump]>0);
ihi = isupp([Stim.ProfileLowFreq]==Lowprofile(4) & [Stim.ProfileSPLjump]>0);
iref = isupp([Stim.ProfileSPLjump]==0);
AAlo = apple(DS(ilo),2);
AAml = apple(DS(iml),2);
AAmh = apple(DS(imh),2);
AAhi = apple(DS(ihi),2);
AAref = apple(DS(iref),2)
AA = [AAlo AAml AAmh AAhi]; % all 1-comp supp data
Fprofile = unique([AA.Fprofile])
baseSPL = unique([AA.baseSPL])
SPLjump = unique([AA.SPLjump])
Fprim=DS(ilo(1)).Stim.Fzwuis/1e3;
%-
apple(DS(ilo),2,1);  % 6.1-kHz supp
apple(DS(ilo),2,2);  % 6.5-kHz supp
apple(DS(ilo),2,3);  % 7.1-kHz supp
apple(DS(ilo),2,4);  % 7.5-kHz supp
apple(DS(ilo),2,5);  % 8.1-kHz supp
apple(DS(ilo),2,6);  % 8.1-kHz supp
%-
apple(DS(iml),2,1);  % 15.3-kHz supp
apple(DS(iml),2,2);  % 15.8-kHz supp
apple(DS(iml),2,3);  % 16.3-kHz supp
apple(DS(iml),2,4);  % 16.9-kHz supp
apple(DS(iml),2,5);  % 17.3-kHz supp
apple(DS(iml),2,6);  % 17.7-kHz supp
%-
apple(DS(imh),2,1);  % 15.3-kHz supp
apple(DS(imh),2,2);  % 15.8-kHz supp
apple(DS(imh),2,3);  % 16.3-kHz supp
apple(DS(imh),2,4);  % 16.9-kHz supp
apple(DS(imh),2,5);  % 17.3-kHz supp
apple(DS(imh),2,6);  % 17.7-kHz supp
%-
apple(DS(ihi),2,1);  % 18.3-kHz supp
apple(DS(ihi),2,2);  % 18.8-kHz supp
apple(DS(ihi),2,3);  % 19.5-kHz supp
apple(DS(ihi),2,4);  % 19.9-kHz supp
apple(DS(ihi),2,5);  % 20.5-kHz supp
apple(DS(ihi),2,6);  % 20.8-kHz supp
%-

%======RG12401===============================
DS = getDS('RG12401')
DX = 200; % um distance between beads 
Lw = CVlocalTRF(DS(36:44), DS(17:25),2, DX); % wideband
L3 = CVlocalTRF(DS(50:54), DS(31:35),2, DX); %=3-7 kHz
L6 = CVlocalTRF(DS(56:60), DS(82:86),2, DX); %=6-10 kHz
L9 = CVlocalTRF(DS(61:65), DS(77:81),2, DX); %=9-13 kHz
L12 = CVlocalTRF(DS(45:49), DS(26:30),2, DX); %=12-16 kHz
L15 = CVlocalTRF(DS(66:70), DS(72:76),2, DX); %=15-19 kHz
save(fullfile(processed_datadir, 'CV\RG12401', 'localTrf_401.mat'), 'DX', 'Lw', 'L3', 'L6', 'L9', 'L12', 'L15');
% phase vs Gain
L  = Lw; Dphase = @(ii)double(ii==3 | ii==9);
L  = L12; Dphase = @(ii)0;
L  = L9; Dphase = @(ii)0;%double(ii==3 | ii==9);
L  = L6; Dphase = @(ii)0;%double(ii==3 | ii==9);
L  = L3; Dphase = @(ii)double(ii<=3);
L  = L15; Dphase = @(ii)-double(ii<=2);
%figure; 
set(gcf,'units', 'normalized', 'position', [0.354 0.314 0.389 0.591]);
CM = jet(numel(L.A2B));
for ii=1:numel(L.A2B); 
    ha(ii) = xplot(-Dphase(ii)-L.A2B(ii).Phase+pmask(L.A2B(ii).Alpha<=0.001), L.A2B(ii).Gain, ...
        'v-', 'color', CM(ii,:), 'markersize', 4, 'linewidth', 3); 
end
grid on
for ddx=0.5:0.2:1,1, xplot(-ylim/a2db(exp(2*pi))+ddx, ylim, 'k-'); end
legend(ha, {L.A2B.LegStr}, 'location', 'northeast', 'fontsize', 10)
set(gca, 'fontsize',12)
ylabel('local BM gain (dB)', 'fontsize', 13); xlabel('Minus local BM phase (cycle)', 'fontsize', 13);
ylim([-22 20]); xlim([0 1.4]);
% =========== note added 19-Mar-2012 14:18:00=========
% Calibration from stapes
% ====================================================%  1     1-1-FS    24k10:44k7 Hz   65 dB        no mod                    85x2x700/1000 ms         L    --     AD-1&2   
%  2     1-2-FS    5k:24k10 Hz     80 dB        no mod                    233x2x700/1000 ms        L    --     AD-1&2   
%  3     1-3-FS    51:5k Hz        80 dB        no mod                    166x2x700/1000 ms ***    L    --     AD-1&2   
%  4     1-4-ZW    0:0 Hz          80 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-1&2   
%  5     1-5-ZW    0:0 Hz          70 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-1&2   
%  6     1-6-ZW    0:0 Hz          60 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-1&2   
%  7     1-7-ZW    0:0 Hz          50 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-1&2   
%  8     1-8-ZW    0:0 Hz          40 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-1&2   
%  9     1-9-ZW    0:0 Hz          30 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-1&2   
%  10    1-10-ZW   0:0 Hz          20 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-1&2   
%  11    1-11-ZW   0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-1&2   
%  12    1-12-ZW   0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-1&2   
% ========Measurements for bead 1
% XXX 13    2-1-ZW    0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
% XXX 14    2-2-ZW    0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
% XXX 15    2-3-ZW    0:0 Hz          20 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
% XXX 16    2-4-ZW    0:0 Hz          30 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=wideband
%  17    3-1-ZW    0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  18    3-2-ZW    0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  19    3-3-ZW    0:0 Hz          20 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  20    3-4-ZW    0:0 Hz          30 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  21    3-5-ZW    0:0 Hz          40 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  22    3-6-ZW    0:0 Hz          50 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  23    3-7-ZW    0:0 Hz          60 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  24    3-8-ZW    0:0 Hz          70 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  25    3-9-ZW    0:0 Hz          80 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=12-16 kHz
%  26    3-10-ZW   0:0 Hz          0 dB/cmp     12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  27    3-11-ZW   0:0 Hz          20 dB/cmp    12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  28    3-12-ZW   0:0 Hz          40 dB/cmp    12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  29    3-13-ZW   0:0 Hz          60 dB/cmp    12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  30    3-14-ZW   0:0 Hz          80 dB/cmp    12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=3-7 kHz
%  31    3-15-ZW   0:0 Hz          0 dB/cmp     3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  32    3-16-ZW   0:0 Hz          20 dB/cmp    3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  33    3-17-ZW   0:0 Hz          40 dB/cmp    3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  34    3-18-ZW   0:0 Hz          80 dB/cmp    3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  35    3-19-ZW   0:0 Hz          60 dB/cmp    3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
% ====================================================
%Picture captured: D:\RawData\EARLY\expdata\corst\RG12401\pictures\bead2.tif
% =========== note added 19-Mar-2012 16:38:28=========
% picture "bead2" shows bead 1 as well: bead 1 is to the lower left; bead 2 is in the middle
%==wideband
%  36    4-1-ZW    0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  37    4-2-ZW    0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  38    4-3-ZW    0:0 Hz          20 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  39    4-4-ZW    0:0 Hz          30 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  40    4-5-ZW    0:0 Hz          40 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  41    4-6-ZW    0:0 Hz          50 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  42    4-7-ZW    0:0 Hz          60 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  43    4-8-ZW    0:0 Hz          70 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  44    4-9-ZW    0:0 Hz          80 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=12-16 kHz
%  45    4-10-ZW   0:0 Hz          0 dB/cmp     12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  46    4-11-ZW   0:0 Hz          20 dB/cmp    12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  47    4-12-ZW   0:0 Hz          40 dB/cmp    12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  48    4-13-ZW   0:0 Hz          60 dB/cmp    12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  49    4-14-ZW   0:0 Hz          80 dB/cmp    12000-16000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=3-7 kHz
%XXX  50    4-15-ZW   0:0 Hz          0 dB/cmp     3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  51    4-16-ZW   0:0 Hz          20 dB/cmp    3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  52    4-17-ZW   0:0 Hz          40 dB/cmp    3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  53    4-18-ZW   0:0 Hz          60 dB/cmp    3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  54    4-19-ZW   0:0 Hz          80 dB/cmp    3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  55    4-20-ZW   0:0 Hz          0 dB/cmp     3000-7000 Hz     30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=6-10 kHz
%  56    4-21-ZW   0:0 Hz          0 dB/cmp     6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  57    4-22-ZW   0:0 Hz          20 dB/cmp    6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  58    4-23-ZW   0:0 Hz          40 dB/cmp    6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  59    4-24-ZW   0:0 Hz          60 dB/cmp    6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  60    4-25-ZW   0:0 Hz          80 dB/cmp    6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=9-13 kHz
%  61    4-26-ZW   0:0 Hz          0 dB/cmp     9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  62    4-27-ZW   0:0 Hz          20 dB/cmp    9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  63    4-28-ZW   0:0 Hz          40 dB/cmp    9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  64    4-29-ZW   0:0 Hz          60 dB/cmp    9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  65    4-30-ZW   0:0 Hz          80 dB/cmp    9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=15-19 kHz
%  66    4-31-ZW   0:0 Hz          0 dB/cmp     15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  67    4-32-ZW   0:0 Hz          20 dB/cmp    15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  68    4-33-ZW   0:0 Hz          40 dB/cmp    15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  69    4-34-ZW   0:0 Hz          60 dB/cmp    15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  70    4-35-ZW   0:0 Hz          80 dB/cmp    15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
% ====================================================
% =========== note added 19-Mar-2012 17:09:11=========
% Back to bead 1
% wideband check
%  71    5-1-ZW    0:0 Hz          20 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=15-19 kHz
%  72    5-2-ZW    0:0 Hz          0 dB/cmp     15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  73    5-3-ZW    0:0 Hz          20 dB/cmp    15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  74    5-4-ZW    0:0 Hz          40 dB/cmp    15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  75    5-5-ZW    0:0 Hz          60 dB/cmp    15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  76    5-6-ZW    0:0 Hz          80 dB/cmp    15000-19000 Hz   30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=9-13 kHz
%  77    5-7-ZW    0:0 Hz          0 dB/cmp     9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  78    5-8-ZW    0:0 Hz          20 dB/cmp    9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  79    5-9-ZW    0:0 Hz          40 dB/cmp    9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  80    5-10-ZW   0:0 Hz          60 dB/cmp    9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  81    5-11-ZW   0:0 Hz          80 dB/cmp    9000-13000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%=6-10 kHz
%  82    5-12-ZW   0:0 Hz          0 dB/cmp     6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  83    5-13-ZW   0:0 Hz          20 dB/cmp    6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  84    5-14-ZW   0:0 Hz          40 dB/cmp    6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  85    5-15-ZW   0:0 Hz          60 dB/cmp    6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  86    5-16-ZW   0:0 Hz          80 dB/cmp    6000-10000 Hz    30 cmp   1x1x15050/15550 ms       L    --     AD-2     
% wideband check
%  87    5-17-ZW   0:0 Hz          20 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
%  88    5-18-ZW   0:0 Hz          20 dB/cmp    100-25000 Hz     40 cmp   1x1x15050/15550 ms       L    --     AD-2     
% ====================================================
% =========== note added 19-Mar-2012 19:44:11=========
% last few recs (89-100) to check on degrading prep.
%  89    5-19-ZW   0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  90    5-20-ZW   0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  91    5-21-ZW   0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  92    5-22-ZW   0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  93    5-23-ZW   0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  94    5-24-ZW   0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  95    5-25-ZW   0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  96    5-26-ZW   0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  97    5-27-ZW   0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  98    5-28-ZW   0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  99    5-29-ZW   0:0 Hz          10 dB/cmp    100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2     
%  100   5-30-ZW   0:0 Hz          0 dB/cmp     100-25000 Hz     40 cmp   1x1x115050/115550.1 ms   L    --     AD-2 Picture captured: D:\RawData\EARLY\expdata\corst\RG12401\pictures\bead1.tif
DS = getDS('RG12402')
apple(DS(6:14),2,0) % wideband
apple(DS(15:23),2,0) % 12-16 kHz
apple(DS(24:31),2,0) % 15-21 kHz







