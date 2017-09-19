function Y = twobeadstuff(kw, varargin);
% twobeadstuff - preliminary processing & overview of 2-bead data
%   twobeadstuff applesave RG12401
%   obtains local transfer functions of experiment RG12401 and saves them
%   the processed_datadir.
%
%   twobeadstuff('getapple', 'RG12401') returns the previously saved 
%   transfer functions.
%
%   twobeadstuff('show', 'RG12401') plots a number of elementary analyses.
%
%   twobeadstuff('CFmatch', 'RG12401', dA, freqRatio) plots single-bead
%   data with shifted amplitude and scaled frequency. This facilitates the
%   estimation of CFs and their ratio.
%
%   To consult Corstiaen's notes:
%   edit \\clust\processed_data$\CV\TwoBeadAnalyses.m
%
%  Quote from these notes
%  The "good" ones are:
%   RG12401 --> perfect; some TTS
%   RG12405 --> low SPL and <5 kHz are strange; possibly some TTS
%   RG12407 --> strange <5 kHz
%   RG12417 --> good, strange <5 kHz
%   RG12422 --> perfect!
%   RG12423 --> good; look at 3rd bead
%   RG12424 --> good
%=========simplified stim========
%   RG12432 --> good
%   RG12433 --> excellent
%   RG12434 --> good
%   RG12436 --> good
%   RG12438 --> three beads
%
%  To adjust slipped cycles in phase:
%   qq = twobeadstuff('getapple', 'RG12424'); aa; CVlocalTRF(qq.Lw, 'plotns'); xlog125([1 35])
%   iend=5; qq = twobeadstuff('getapple', 'RG12424'); aa; for jj=1:iend, L=qq.(qq.LnNames{jj}); L.A2B = L.A2B(1:end); CVlocalTRF(L, 'plotns'); end; xlog125([1 35])


switch lower(kw),
    case 'saveall',
        EXP = {'RG12401' 'RG12405' 'RG12407' 'RG12417' 'RG12422' 'RG12423' 'RG12424' ...
            'RG12432' 'RG12433' 'RG12434' 'RG12436' };
        for ii=1:numel(EXP), local_save(EXP{ii}); end;
    case 'applesave',
        local_save(varargin{1});
    case 'getapple',
        Y = local_get(varargin{1});
    case 'show',
        local_show(varargin{1});
    case 'cfmatch',
        local_cfmatch(varargin{:});
end



%======================

function Sdir = local_Sdir();
Sdir = fullfile(processed_datadir, 'CV', 'twobeads');


function local_save(ExpName);
Sfile = fullfile(local_Sdir, [ExpName '_localTRF.mat']);
noPlot = false;
DS = getds(ExpName);
switch upper(ExpName),
    case 'RG12401',
        ME = ds2trf(DS(1:3),2,0); % calibration from stapes
        AC = ds2trf(DS(1:3),1,0); % calibration from mike
        % Distance between beads; radial position is very similar for both beads
        dX = diff([220 400]);
        dY = diff([400 365]);
        dZ = diff([0 0]);
        BeadDist = sqrt(sum([dX dY dZ].^2));
        % BEAD 1
        % CF = 13.42
        % Indices to groups of recs
        B1iBB = 17:25; % 0.1-25 kHz
        B1i15_19 = 72:76; % 15k-19k
        B1i12_16 = 26:30; % 12k-16k
        B1i9_13 = 77:81; % 9k-13k
        B1i6_10 = 82:86; % 6k-10k
        B1i3_7 = 31:35; % 3k-7k
        B1prephealth = [71 87:100]; % 0.1-25 kHz to check on prep health
        % apples
        AP_Aw = apple(DS(B1iBB),2,0,ME)
        AP_A15 = apple(DS(B1i15_19),2,0,ME)
        AP_A12 = apple(DS(B1i12_16),2,0,ME)
        AP_A9 = apple(DS(B1i9_13),2,0,ME)
        AP_A6 = apple(DS(B1i6_10),2,0,ME)
        AP_A3 = apple(DS(B1i3_7),2,0,ME)
        AP_Ah = apple(DS(B1prephealth),2,0,ME)
        % BEAD 2
        % CF = 14.75
        % Indices to groups of recs
        B2iBB = 36:44; % 0.1-25 kHz
        B2i15_19 = 66:70; % 15k-19k
        B2i12_16 = 45:49; % 12k-16k
        B2i9_13 = 61:65; % 9k-13k
        B2i6_10 = 56:60; % 6k-10k
        B2i3_7 = 51:55; % 3k-7k
        % apples
        AP_Bw = apple(DS(B2iBB),2,0,ME)
        AP_B15 = apple(DS(B2i15_19),2,0,ME)
        AP_B12 = apple(DS(B2i12_16),2,0,ME)
        AP_B9 = apple(DS(B2i9_13),2,0,ME)
        AP_B6 = apple(DS(B2i6_10),2,0,ME)
        AP_B3 = apple(DS(B2i3_7),2,0,ME)
        % TRANSFER FUNCTIONS
        Lw = CVlocalTRF(DS(B2iBB),DS(B1iBB),2,BeadDist,noPlot) % wideband
        L15 = CVlocalTRF(DS(B2i15_19),DS(B1i15_19),2,BeadDist,noPlot) % 15k-19k
        L12 = CVlocalTRF(DS(B2i12_16),DS(B1i12_16),2,BeadDist,noPlot) % 12k-16k
        L9 = CVlocalTRF(DS(B2i9_13),DS(B1i9_13),2,BeadDist,noPlot) % 9k-13k
        L6 = CVlocalTRF(DS(B2i6_10),DS(B1i6_10),2,BeadDist,noPlot) % 6k-10k
        L3 = CVlocalTRF(DS(B2i3_7),DS(B1i3_7),2,BeadDist,noPlot) % 3k-7k
        LnNames = {'L3' 'L6' 'L9' 'L12' 'L15' };
    case 'RG12405',
        ME = ds2trf(DS(1:3),2,0); % calibration from stapes
        AC = ds2trf(DS(1:3),1,0); % calibration from mike
        % Distance between beads
        dX = diff([84 283]);
        dY = diff([388 458]);
        dZ = diff([37175.1 36882.6]);
        BeadDist = sqrt(sum([dX dY dZ].^2));
        % BEAD 1
        % CF = 12.17
        % Indices to groups of recs
        B1iBB = [4:7:32 109:7:130]; % 0.1-25 kHz
        B1i15_19 = [5:7:33 110:7:131]; % 15k-19k
        B1i12_16 = [6:7:34 111:7:132]; % 12k-16k
        B1i9_13 = [7:7:35 112:7:133]; % 9k-13k
        B1i6_10 = [8:7:36 113:7:134]; % 6k-10k
        B1i3_7 = [9:7:37 114:7:135]; % 3k-7k
        B1i1_4 = [10:7:38 115:7:136]; % 0.1k-4k
        B1prephealth = [4 18 32 39:41 107:108 137:142]; % 0.1-25 kHz to check on prep health
        % BEAD 2
        % CF = 14.21
        % Indices to groups of recs
        B2iBB = 42:7:98; % 0.1-25 kHz
        B2i15_19 = 43:7:99; % 15k-19k
        B2i12_16 = 44:7:100; % 12k-16k
        B2i9_13 = 45:7:101; % 9k-13k
        B2i6_10 = 46:7:102; % 6k-10k
        B2i3_7 = 47:7:103; % 3k-7k
        B2i1_4 = 48:7:104; % 0.1k-4k
        B2prephealth = [56 105:106]; % 0.1-25 kHz to check on prep health
        % TRANSFER FUNCTIONS: does DX have to be longitudinal distance only..?
        Lw = CVlocalTRF(DS(B2iBB),DS(B1iBB),2,BeadDist,noPlot) % wideband
        L15 = CVlocalTRF(DS(B2i15_19),DS(B1i15_19),2,BeadDist,noPlot) % 15k-19k
        L12 = CVlocalTRF(DS(B2i12_16),DS(B1i12_16),2,BeadDist,noPlot) % 12k-16k
        L9 = CVlocalTRF(DS(B2i9_13),DS(B1i9_13),2,BeadDist,noPlot) % 9k-13k
        L6 = CVlocalTRF(DS(B2i6_10),DS(B1i6_10),2,BeadDist,noPlot) % 6k-10k
        L3 = CVlocalTRF(DS(B2i3_7),DS(B1i3_7),2,BeadDist,noPlot) % 3k-7k
        L1 = CVlocalTRF(DS(B2i1_4),DS(B1i1_4),2,BeadDist,noPlot) % 0.1k-4k
        LnNames = {'L1' 'L3' 'L6' 'L9' 'L12' 'L15' };
    case 'RG12407',
        % Initialization
        ME = ds2trf(DS(1:3),2,0); % calibration from stapes
        AC = ds2trf(DS(1:3),1,0); % calibration from mike
        % Distance between beads
        dX = diff([370 410]);
        dY = diff([64 170]);
        dZ = diff([23983.5 23853.5]);
        BeadDist = sqrt(sum([dX dY dZ].^2));
        % BEAD 1
        % Single bead
        % CF = 11.64
        % Indices to groups of recs
        B1iBB = [4 25 32 53 60 81 88 109 116]; % 0.1-25 kHz
        B1i15_19 = [5 26 33 54 61 82 89 110 117]; % 15k-19k
        B1i12_16 = [6 27 34 55 62 83 90 111 118]; % 12k-16k
        B1i9_13 = [7 28 35 56 63 84 91 112 120]; % 9k-13k
        B1i6_10 = [8 29 36 57 64 85 92 113 121]; % 6k-10k
        B1i3_7 = [9 30 37 58 65 86 93 114 122]; % 3k-7k
        B1i1_4 = [10 31 38 59 66 87 94 115 123]; % 0.1k-4k
        % BEAD 2
        % Last bead of row of three
        % CF = 13.42
        % Indices to groups of recs
        B2iBB = [11 18 39 46 67 74 95 102 124]; % 0.1-25 kHz
        B2i15_19 = [12 19 40 47 68 75 96 103 125]; % 15k-19k
        B2i12_16 = [13 20 41 48 69 76 97 104 126]; % 12k-16k
        B2i9_13 = [14 21 42 49 70 77 98 105 127]; % 9k-13k
        B2i6_10 = [15 22 43 50 71 78 99 106 128]; % 6k-10k
        B2i3_7 = [16 23 44 51 72 79 100 107 129]; % 3k-7k
        B2i1_4 = [17 24 45 52 73 80 101 108 130]; % 0.1k-4k
        % TRANSFER FUNCTIONS: does DX have to be longitudinal distance only..?
        Lw = CVlocalTRF(DS(B2iBB),DS(B1iBB),2,BeadDist,noPlot) % wideband
        L15 = CVlocalTRF(DS(B2i15_19),DS(B1i15_19),2,BeadDist,noPlot) % 15k-19k
        L12 = CVlocalTRF(DS(B2i12_16),DS(B1i12_16),2,BeadDist,noPlot) % 12k-16k
        L9 = CVlocalTRF(DS(B2i9_13),DS(B1i9_13),2,BeadDist,noPlot) % 9k-13k
        L6 = CVlocalTRF(DS(B2i6_10),DS(B1i6_10),2,BeadDist,noPlot) % 6k-10k
        L3 = CVlocalTRF(DS(B2i3_7),DS(B1i3_7),2,BeadDist,noPlot) % 3k-7k
        L1 = CVlocalTRF(DS(B2i1_4),DS(B1i1_4),2,BeadDist,noPlot) % 0.1k-4k
        LnNames = {'L1' 'L3' 'L6' 'L9' 'L12' 'L15' };
    case 'RG12417',
        % 10x objective
        % Though this prep should hve been stable (see also pictures), bead 1 does
        % show strange behavior in LF tail, similar to a drifting bead.
        % Initialization
        ME = ds2trf(DS(185:186),2,0); % calibration from stapes
        AC = ds2trf(DS(185:186),1,0); % calibration from mike
        % Distance between beads 2 and 3; A LARGE DIFFERENCE IN RADIAL POSITION
        % EXISTS IN THIS CASE!!
        dX = diff([48 42]);
        dY = diff([419 253]);
        dZ = diff([0 -861.5]);
        BeadDist = sqrt(sum([dX dY dZ].^2));
        % BEAD 1
        % CF = 16.7
        % Bead was discontinued after 40 dB due to drift or group of beads being in
        % focus and not returning to same location each time
        % Indices to groups of recs
        B1iBB = [1 22 29 50 57]; % 0.1-25 kHz
        B1i15_19 = [2 23 30 51 58]; % 15k-19k
        B1i12_16 = [3 24 31 52 59]; % 12k-16k
        B1i9_13 = [4 25 32 53 60]; % 9k-13k
        B1i6_10 = [5 26 33 54 61]; % 6k-10k
        B1i3_7 = [6 27 34 55 62]; % 3k-7k
        B1i1_4 = [7 28 35 56 63]; % 0.1k-4k
        % BEAD 2
        % Individual bead
        % CF = 17.22
        % Indices to groups of recs
        B2iBB = [8 15 36 43 64 71 120 141 148]; % 0.1-25 kHz
        B2i15_19 = [9 16 37 44 65 72 121 142 149]; % 15k-19k
        B2i12_16 = [10 17 38 45 66 73 122 143 150]; % 12k-16k
        B2i9_13 = [11 18 39 46 67 74 123 144 151]; % 9k-13k
        B2i6_10 = [12 19 40 47 68 75 124 145 152]; % 6k-10k
        B2i3_7 = [13 20 41 48 69 76 125 146 153]; % 3k-7k
        B2i1_4 = [14 21 42 49 70 77 126 147 154]; % 0.1k-4k
        % BEAD 3
        % Individual bead
        % CF = 16.12
        % Indices to groups of recs
        B3iBB = [78 85 92 99 106 113 127 134 155]; % 0.1-25 kHz
        B3i15_19 = [79 86 93 100 107 114 128 135 156]; % 15k-19k
        B3i12_16 = [80 87 94 101 108 115 129 136 157]; % 12k-16k
        B3i9_13 = [81 88 95 102 109 116 130 137 158]; % 9k-13k
        B3i6_10 = [82 89 96 103 110 117 131 138 159]; % 6k-10k
        B3i3_7 = [83 90 97 104 111 118 132 139 160]; % 3k-7k
        B3i1_4 = [84 91 98 105 112 119 133 140 161]; % 0.1k-4k
        % TRANSFER FUNCTIONS: does DX have to be longitudinal distance only..?
        Lw = CVlocalTRF(DS(B2iBB),DS(B3iBB),2,BeadDist,noPlot) % wideband
        L15 = CVlocalTRF(DS(B2i15_19),DS(B3i15_19),2,BeadDist,noPlot) % 15k-19k
        L12 = CVlocalTRF(DS(B2i12_16),DS(B3i12_16),2,BeadDist,noPlot) % 12k-16k
        L9 = CVlocalTRF(DS(B2i9_13),DS(B3i9_13),2,BeadDist,noPlot) % 9k-13k
        L6 = CVlocalTRF(DS(B2i6_10),DS(B3i6_10),2,BeadDist,noPlot) % 6k-10k
        L3 = CVlocalTRF(DS(B2i3_7),DS(B3i3_7),2,BeadDist,noPlot) % 3k-7k
        L1 = CVlocalTRF(DS(B2i1_4),DS(B3i1_4),2,BeadDist,noPlot) % 0.1k-4k
        LnNames = {'L1' 'L3' 'L6' 'L9' 'L12' 'L15' };
    case 'RG12422',
        % 10x objective
        % For this prep special care was taken to have two solitary beads before
        % starting the measurements. This succeeded and resulted in a very good
        % experiment.
        % Initialization
        ME = ds2trf(DS(144:146),2,0); % calibration from stapes
        AC = ds2trf(DS(144:146),1,0); % calibration from mike
        % Distance between beads; bead 1 closer to OSL than bead 2 (ca. 1/5 vs 1/2 of BM)
        dX = diff([498 354]);
        dY = diff([319 185]);
        dZ = diff([22428.7 22573.1]);
        BeadDist = sqrt(sum([dX dY dZ].^2));
        % BEAD 1
        % CF = 21.68
        % Indices to groups of recs
        B1iBB = [2 20 26 44 50 68 74 92 98]; % 0.1-30 kHz
        B1i22_30 = [3 21 27 45 51 69 75 93]; % 22k-30k
        B1i16_24 = [4 22 28 46 52 70 76 94 99]; % 16k-24k
        B1i10_18 = [5 23 29 47 53 71 77 95 100]; % 10k-18k
        B1i4_12 = [6 24 30 48 54 72 78 96 101]; % 4k-12k
        B1i1_8 = [7 25 31 49 55 73 79 97 102]; % 100-8k
        % apples
        AP_Aw = apple(DS(B1iBB),2,0,ME)
        AP_A22 = apple(DS(B1i22_30),2,0,ME)
        AP_A16 = apple(DS(B1i16_24),2,0,ME)
        AP_A10 = apple(DS(B1i10_18),2,0,ME)
        AP_A4 = apple(DS(B1i4_12),2,0,ME)
        AP_A1 = apple(DS(B1i1_8),2,0,ME)
        % BEAD 2
        % Individual bead
        % CF = 18.63
        % Indices to groups of recs
        B2iBB = [8 14 32 38 56 61 80 86 103]; % 0.1-30 kHz
        B2i22_30 = [9 15 33 39 57 63 81 87]; % 22k-30k
        B2i16_24 = [10 16 34 40 58 64 82 88 104]; % 16k-24k
        B2i10_18 = [11 17 35 41 59 65 83 89 105]; % 10k-18k
        B2i4_12 = [12 18 36 42 60 66 84 90 106]; % 4k-12k
        B2i1_8 = [13 19 37 43 62 67 85 91 107]; % 100-8k
        % apples
        AP_Bw = apple(DS(B2iBB),2,0,ME)
        AP_B22 = apple(DS(B2i22_30),2,0,ME)
        AP_B16 = apple(DS(B2i16_24),2,0,ME)
        AP_B10 = apple(DS(B2i10_18),2,0,ME)
        AP_B4 = apple(DS(B2i4_12),2,0,ME)
        AP_B1 = apple(DS(B2i1_8),2,0,ME)
        % TRANSFER FUNCTIONS: does DX have to be longitudinal distance only..?
        Lw = CVlocalTRF(DS(B1iBB),DS(B2iBB),2,BeadDist,noPlot) % wideband
        L22 = CVlocalTRF(DS(B1i22_30),DS(B2i22_30),2,BeadDist,noPlot) % 22k-30k
        L16 = CVlocalTRF(DS(B1i16_24),DS(B2i16_24),2,BeadDist,noPlot) % 16k-24k
        L10 = CVlocalTRF(DS(B1i10_18),DS(B2i10_18),2,BeadDist,noPlot) % 10k-18k
        L4 = CVlocalTRF(DS(B1i4_12),DS(B2i4_12),2,BeadDist,noPlot) % 4k-12k
        L1 = CVlocalTRF(DS(B1i1_8),DS(B2i1_8),2,BeadDist,noPlot) % 1k-8k
        LnNames = {'L1' 'L4' 'L10' 'L16' 'L22'};
    case 'RG12423',
        %% RG12423
        % 5x objective
        % This prep was tried as a three-bead prep. However, reflectance of beads 1
        % and 2 was very poor, causing the experiment to take a lot of time and
        % possibly have health degraded. Bead 2 seems to have drifted. Also, the 70 and 80 dB recs were not
        % completed due to loss of beads 1 and 2.
        % Initialization
        ME = ds2trf(DS(147:149),2,0); % calibration from stapes
        AC = ds2trf(DS(147:149),1,0); % calibration from mike
        % Distance between beads; bead 1 ca. 1/4 of BM; bead 2 ca. 1/5 of BM; bead 3 ca. 1/4 of BM
        dX1_2 = diff([273 234]);
        dY1_2 = diff([354 313]);
        dZ1_2 = diff([25448.1 25555.3]);
        dist1_2 = sqrt(sum([dX1_2 dY1_2 dZ1_2].^2));
        dX2_3 = diff([234 158]);
        dY2_3 = diff([313 255]);
        dZ2_3 = diff([25555.3 25613.8]);
        dist2_3 = sqrt(sum([dX2_3 dY2_3 dZ2_3].^2));
        dX1_3 = diff([273 158]);
        dY1_3 = diff([354 255]);
        dZ1_3 = diff([25448.1 25613.8]);
        dist1_3 = sqrt(sum([dX1_3 dY1_3 dZ1_3].^2));
        % BEAD 1
        % CF = 18.63
        % Individual bead
        % Indices to groups of recs
        B1iBB = [12 44 50 80 86 116 122]; % 0.1-30 kHz
        B1i22_30 = [13 45 51 81 87 117 123]; % 22k-30k
        B1i16_24 = [14 46 52 82 88 119 124]; % 16k-24k
        B1i10_18 = [15 47 53 83 89 118 125]; % 10k-18k
        B1i4_12 = [16 48 54 84 90 120 126]; % 4k-12k
        B1i1_8 = [17 49 55 85 91 121 127]; % 100-8k
        % apples
        AP_Aw = apple(DS(B1iBB),2,0,ME)
        AP_A22 = apple(DS(B1i22_30),2,0,ME)
        AP_A16 = apple(DS(B1i16_24),2,0,ME)
        AP_A10 = apple(DS(B1i10_18),2,0,ME)
        AP_A4 = apple(DS(B1i4_12),2,0,ME)
        AP_A1 = apple(DS(B1i1_8),2,0,ME)
        % BEAD 2
        % Bead has drifted!
        % Individual bead
        % CF = 16.82
        % Indices to groups of recs
        B2iBB = [19 37 43 56 74 92 110 128]; % 0.1-30 kHz
        B2i22_30 = [20 38 57 75 93 111 129]; % 22k-30k
        B2i16_24 = [21 39 58 76 94 112 130]; % 16k-24k
        B2i10_18 = [22 40 59 77 95 113 131]; % 10k-18k
        B2i4_12 = [23 41 60 78 96 114 132]; % 4k-12k
        B2i1_8 = [24 42 61 79 97 115 133]; % 100-8k
        % BEAD 3
        % Individual bead
        % CF = 16.16
        % Indices to groups of recs
        B3iBB = [25 31 62 68 98 104 134 140]; % 0.1-30 kHz
        B3i22_30 = [26 32 63 69 99 105 135 141]; % 22k-30k
        B3i16_24 = [27 33 64 70 100 106 136 142]; % 16k-24k
        B3i10_18 = [28 34 65 71 101 107 137 143]; % 10k-18k
        B3i4_12 = [29 35 66 72 102 108 138 144]; % 4k-12k
        B3i1_8 = [30 36 67 73 103 109 139 145]; % 100-8k
        % apples
        AP_Bw = apple(DS(B3iBB),2,0,ME)
        AP_B22 = apple(DS(B3i22_30),2,0,ME)
        AP_B16 = apple(DS(B3i16_24),2,0,ME)
        AP_B10 = apple(DS(B3i10_18),2,0,ME)
        AP_B4 = apple(DS(B3i4_12),2,0,ME)
        AP_B1 = apple(DS(B3i1_8),2,0,ME)
        %
        % picking bead 1&3; #2 has drifted.
        BeadDist = dist1_3;
        [dX, dY, dZ] = deal(dX1_3, dY1_3, dZ1_3);
        % TRANSFER FUNCTIONS: does DX have to be longitudinal distance only..?
        %         figure, CVlocalTRF(DS(B1iBB),DS(B2iBB),2,dist1_2,1) % wideband
        %         figure, CVlocalTRF(DS(B1i22_30),DS(B2i22_30),2,dist1_2,1) % 22k-30k
        %         figure, CVlocalTRF(DS(B1i16_24),DS(B2i16_24),2,dist1_2,1) % 16k-24k
        %         figure, CVlocalTRF(DS(B1i10_18),DS(B2i10_18),2,dist1_2,1) % 10k-18k
        %         figure, CVlocalTRF(DS(B1i4_12),DS(B2i4_12),2,dist1_2,1) % 4k-12k
        %         figure, CVlocalTRF(DS(B1i1_8),DS(B2i1_8),2,dist1_2,1) % 1k-8k
        %         figure, CVlocalTRF(DS([B1i1_8 B1i4_12 B1i10_18 B1i16_24 B1i22_30]),DS([B2i1_8 B2i4_12 B2i10_18 B2i16_24 B2i22_30]),2,dist1_2,1) % all
        %         %
        %         figure, CVlocalTRF(DS(B2iBB),DS(B3iBB),2,dist2_3,1) % wideband
        %         figure, CVlocalTRF(DS(B2i22_30),DS(B3i22_30),2,dist2_3,1) % 22k-30k
        %         figure, CVlocalTRF(DS(B2i16_24),DS(B3i16_24),2,dist2_3,1) % 16k-24k
        %         figure, CVlocalTRF(DS(B2i10_18),DS(B3i10_18),2,dist2_3,1) % 10k-18k
        %         figure, CVlocalTRF(DS(B2i4_12),DS(B3i4_12),2,dist2_3,1) % 4k-12k
        %         figure, CVlocalTRF(DS(B2i1_8),DS(B3i1_8),2,dist2_3,1) % 1k-8k
        %         figure, CVlocalTRF(DS([B2i1_8 B2i4_12 B2i10_18 B2i16_24 B2i22_30]),DS([B3i1_8 B3i4_12 B3i10_18 B3i16_24 B3i22_30]),2,dist2_3,1) % all
        %
        Lw = CVlocalTRF(DS(B1iBB),DS(B3iBB),2,BeadDist,noPlot) % wideband
        L22 = CVlocalTRF(DS(B1i22_30),DS(B3i22_30),2,BeadDist,noPlot) % 22k-30k
        L16 = CVlocalTRF(DS(B1i16_24),DS(B3i16_24),2,BeadDist,noPlot) % 16k-24k
        L10 = CVlocalTRF(DS(B1i10_18),DS(B3i10_18),2,BeadDist,noPlot) % 10k-18k
        L4 = CVlocalTRF(DS(B1i4_12),DS(B3i4_12),2,BeadDist,noPlot) % 4k-12k
        L1 = CVlocalTRF(DS(B1i1_8),DS(B3i1_8),2,BeadDist,noPlot) % 1k-8k
        LnNames = {'L1' 'L4' 'L10' 'L16' 'L22'};
    case 'RG12424',
        % 5x objective
        % Both beads were solitary and had same radial position; ca. 1/7 of BM.
        % After tearing the RW, a lot of fluid streamed from cochlea. Though there
        % was no apparent drift on the monitor, the traces of bead 1 did not align
        % the way they should. At the end of the experiment bead 1 drifted a lot
        % all of a sudden.
        % Initialization
        ME = ds2trf(DS(140:142),2,0); % calibration from stapes
        AC = ds2trf(DS(140:142),1,0); % calibration from mike
        % Distance between beads; bead 1 ca. 1/4 of BM; bead 2 ca. 1/5 of BM; bead 3 ca. 1/4 of BM
        dX = diff([358 215]);
        dY = diff([550 413]);
        dZ = diff([29837.8 29930.6]);
        BeadDist = sqrt(sum([dX dY dZ].^2));
        %==== BEAD 1
        % CF = 20.73
        % Bead might have drifted
        % Individual bead
        % Indices to groups of recs
        B1iBB = [1 19 25 43 49 67 73 91 97]; % 0.1-30 kHz
        B1i22_30 = [2 20 26 44 50 68 74 92]; % 22k-30k
        B1i16_24 = [3 21 27 45 51 69 75 93 98]; % 16k-24k
        B1i10_18 = [4 22 28 46 52 70 76 94 99]; % 10k-18k
        B1i4_12 = [5 23 29 47 53 71 77 95 100]; % 4k-12k
        B1i1_8 = [6 24 30 48 54 72 78 96 101]; % 100-8k
        % apples
        AP_Aw = apple(DS(B1iBB),2,0,ME);
        AP_A22 = apple(DS(B1i22_30),2,0,ME);
        AP_A16 = apple(DS(B1i16_24),2,0,ME);
        AP_A10 = apple(DS(B1i10_18),2,0,ME);
        AP_A4 = apple(DS(B1i4_12),2,0,ME);
        AP_A1 = apple(DS(B1i1_8),2,0,ME);
        %===== BEAD 2
        % Individual bead
        % CF = 17.06
        % Indices to groups of recs
        B2iBB = [7 13 31 37 55 61 79 85 102]; % 0.1-30 kHz
        B2i22_30 = [8 14 32 38 56 62 80 86]; % 22k-30k
        B2i16_24 = [9 15 33 39 57 63 81 87 103]; % 16k-24k
        B2i10_18 = [10 16 34 40 58 64 82 88 104]; % 10k-18k
        B2i4_12 = [11 17 35 41 59 65 83 89 105]; % 4k-12k
        B2i1_8 = [12 18 36 42 60 66 84 90 106]; % 100-8k
        % apples
        AP_Bw = apple(DS(B2iBB),2,0,ME)
        AP_B22 = apple(DS(B2i22_30),2,0,ME)
        AP_B16 = apple(DS(B2i16_24),2,0,ME)
        AP_B10 = apple(DS(B2i10_18),2,0,ME)
        AP_B4 = apple(DS(B2i4_12),2,0,ME)
        AP_B1 = apple(DS(B2i1_8),2,0,ME)
        % HF plateau after bead 1 had drifted
        B1iHFP = 131:134; % HF plateau bead 1
        B2iHFP = 135:138; % HF plateau bead 2
        % TRANSFER FUNCTIONS: does DX have to be longitudinal distance only..?
        Lw = CVlocalTRF(DS(B1iBB),DS(B2iBB),2,BeadDist,noPlot) % wideband
        L22 = CVlocalTRF(DS(B1i22_30),DS(B2i22_30),2,BeadDist,noPlot) % 22k-30k
        L16 = CVlocalTRF(DS(B1i16_24),DS(B2i16_24),2,BeadDist,noPlot) % 16k-24k
        L10 = CVlocalTRF(DS(B1i10_18),DS(B2i10_18),2,BeadDist,noPlot) % 10k-18k
        L4 = CVlocalTRF(DS(B1i4_12),DS(B2i4_12),2,BeadDist,noPlot) % 4k-12k
        L1 = CVlocalTRF(DS(B1i1_8),DS(B2i1_8),2,BeadDist,noPlot) % 1k-8k
        LHFP = CVlocalTRF(DS(B1iHFP),DS(B2iHFP),2,BeadDist,noPlot) % HF plateau
        LnNames = {'L1' 'L4' 'L10' 'L16' 'L22' 'LHFP'};
    case 'RG12432',
        %% RG12432
        % 5x objective
        % Pretty large spacing
        % Experiment aimed at passive gain study
        % This experiment we realised that we should not switch beads for each SPL,
        % but just record from one and then the other. Data in iRec: [6:19 26:34 43:48 55:59]
        % Initialization
        ME = ds2trf(DS(81:83),2,0); % calibration from stapes
        AC = ds2trf(DS(81:83),1,0); % calibration from mike
        % Distance between beads; bead 1 closer to OSL than bead 2 (ca. 1/5 vs 1/3 of BM)
        dX = diff([40 203]);
        dY = diff([325 92]);
        dZ = diff([27653.3 27896.7]);
        BeadDist = sqrt(sum([dX dY 0*dZ].^2));
        % === BEAD 1
        % CF = ??????
        % Solitary bead); 1/5 of BM
        % Indices to groups of recs
        B1iBB = 1:5; % 0.1-30 kHz
        B1i7_17 = {[22:23 52] [24 53] [35:38 54] 39:42 49:51 79}; % 7k-17k
        % apples
        AP_Aw = apple(DS(B1iBB),2,0,ME)
        AP_An = apple(DS(B1i7_17),2,0,ME)
        % === BEAD 2
        % Double beads; 1/3 of BM
        % CF = ??????
        % Indices to groups of recs
        B2iBB = 69:73; % 0.1-30 kHz
        B2i7_17 = {[64 65 67] [63 66 68] [62 74] [61 75] [60 76] 78}; % 7k-17k
        % apples
        AP_Bw = apple(DS(B2iBB),2,0,ME)
        AP_Bn = apple(DS(B2i7_17),2,0,ME)
        % TRANSFER FUNCTIONS
        Lw = CVlocalTRF(DS(B1iBB),DS(B2iBB),2,BeadDist,noPlot) % wideband
        Ln = CVlocalTRF(DS(B1i7_17),DS(B2i7_17),2,BeadDist,noPlot) % 7k-17k
    case 'RG12433',
        % 5x objective
        % Experiment aimed at passive gain study
        % Initialization
        ME = ds2trf(DS(39:41),2,0); % calibration from stapes
        AC = ds2trf(DS(39:41),1,0); % calibration from mike
        % Distance between beads; beads almost same radial position (1/2 and 3/5 of BM)
        dX = diff([198 84]);
        dY = diff([106 16]);
        dZ = diff([31828.6 31970.2]);
        BeadDist = sqrt(sum([dX dY 0*dZ].^2));
        % ======BEAD 1
        % Changed during experiment; only most recent datataken into account
        % CF = ??????
        % Almost solitary bead; 2 nearby; 1/2 of BM
        % Indices to groups of recs
        B1iBB = {21:22 23 24 25 26 27 29 31}; % 0.1-30 kHz
        B1i8_18 = [16:20 28 30 32]; % 8k-18k
        % apples
        AP_Aw = apple(DS(B1iBB),2,0,ME)
        AP_An = apple(DS(B1i8_18),2,0,ME)
        % ======BEAD 2
        % CF = ??????
        % Almost solitary bead; 2 nearby; 3/5 of BM
        % Indices to groups of recs
        B2iBB = [6:10 33 35 37]; % 0.1-30 kHz
        B2i8_18 = [11:15 34 36 38]; % 8k-18k
        % apples
        AP_Bw = apple(DS(B2iBB),2,0,ME)
        AP_Bn = apple(DS(B2i8_18),2,0,ME)
        % TRANSFER FUNCTIONS
        Lw = CVlocalTRF(DS(B1iBB),DS(B2iBB),2,BeadDist,noPlot) % wideband
        Ln = CVlocalTRF(DS(B1i8_18),DS(B2i8_18),2,BeadDist,noPlot) % 8k-18k
        %=====================================
    case 'RG12434',
        %=====================================
        % 5x objective
        % Experiment aimed at passive gain study
        % Prep possibly had damaged semi-circular canal, no visible effect on data
        % though..
        ME = ds2trf(DS(62:64),2,0); % calibration from stapes
        AC = ds2trf(DS(62:64),1,0); % calibration from mike
        % Distance between beads; beads same radial position (1/3 of BM)
        dX = diff([350 102]);
        dY = diff([173 -41]);
        dZ = diff([26007.6 26122.2]);
        BeadDist = sqrt(sum([dX dY 0*dZ].^2));
        % =====BEAD 1
        % CF = ??????
        % Three bead); 1/3 of BM
        % Indices to groups of recs
        B1iBB = {[1 33:34] [2 30:32] [3 29] [4 28] [5 27] 50 52 54}; % 0.1-30 kHz
        B1i8_18 = {35:37 38:40 42:43 44:46 47:49 51 53 55}; % 8k-18k
        % apples
        AP_Aw = apple(DS(B1iBB),2,0,ME)
        AP_An = apple(DS(B1i8_18),2,0,ME)
        % =====BEAD 2
        % CF = ??????
        % Double bead; 1/3 of BM
        % Indices to groups of recs
        B2iBB = [6:10 56 58 60]; % 0.1-30 kHz
        B2i8_18 = {13 16 17:19 20:22 23:26 57 59 61}; % 8k-18k
        % apples
        AP_Bw = apple(DS(B2iBB),2,0,ME)
        AP_Bn = apple(DS(B2i8_18),2,0,ME)
        % TRANSFER FUNCTIONS
        Lw = CVlocalTRF(DS(B1iBB),DS(B2iBB),2,BeadDist,noPlot) % wideband
        Ln = CVlocalTRF(DS(B1i8_18),DS(B2i8_18),2,BeadDist,noPlot) % 8k-18k
        %=====================================
    case 'RG12436',
        %=====================================
        % 5x objective
        % Experiment aimed at passive gain study
        % Beads close to OSL
        % Initialization
        ME = ds2trf(DS(81:82),2,0); % calibration from stapes; accidentally discarded calibration <5 kHz... (fix that bug, Marcel ;) )
        AC = ds2trf(DS(81:82),1,0); % calibration from mike
        % Distance between beads; beads almost same radial position (1/5 and 1/6 of BM)
        dX = diff([215 7]);
        dY = diff([56 338]);
        dZ = diff([30006.7 30191.8]);
        BeadDist = sqrt(sum([dX dY dZ].^2));
        % BEAD 1
        % CF = ??????
        % Solitary bead; 1/5 of BM; poor reflectivity
        % Indices to groups of recs
        B1iBB = {[4 40:43] [5 44:46] [3 47:49] [2 50:51] [1 52:53] 69 71 73}; % 0.1-30 kHz
        B1i7_17 = {[54:56 60] [57:59 61] 62:64 65:66 67:68 70 72 74}; % 7k5-17k5
        % apples
        AP_Aw = apple(DS(B1iBB),2,0,ME)
        AP_An = apple(DS(B1i7_17),2,0,ME)
        % BEAD 2
        % CF = ??????
        % Solitary bead); 1/6 of BM; seemingly very far into cochlea so angle of laser is "different" (though CF tells differently)
        % Indices to groups of recs
        B2iBB = {38:39 36:37 35 76 75 77 79}; % 0.1-30 kHz
        B2i7_17 = {[20 22 29] [24:25 33] 26:28 30:31 [32 34] 76 78 80}; % 7k5-17k5
        % apples
        AP_Bw = apple(DS(B2iBB),2,0,ME)
        AP_Bn = apple(DS(B2i7_17),2,0,ME)
        % TRANSFER FUNCTIONS
        Lw = CVlocalTRF(DS(B1iBB),DS(B2iBB),2,BeadDist,noPlot) % wideband
        Ln = CVlocalTRF(DS(B1i7_17),DS(B2i7_17),2,BeadDist,noPlot) % 7k5-17k5
        %=====================================
    case {'RG12438_12', 'RG12438_23' 'RG12438_13'},
        %=====================================
        % 5x objective
        % Experiment aimed at passive gain study
        % pretty small spacing
        % accidentally punctured tympanic membrane
        % radial position almost equal for all three beads
        % distance between beads ca. 2*70 um
        %% RG12438 beads 1&2
        % Three bead prep!
        DS =getds('RG12438');
        ME = ds2trf(DS([111 112 114]),2,0); % calibration from stapes
        AC = ds2trf(DS([111 112 114]),1,0); % calibration from mike
        % Distance between beads 1 and 2
        dX12 = diff([309 269]);
        dY12 = diff([80 19]);
        dZ12 = diff([31440.7 31440.7]);
        dist12 = sqrt(sum([dX12 dY12 dZ12].^2));
        % Distance between beads 2 and 3
        dX23 = diff([269 221]);
        dY23 = diff([19 -22]);
        dZ23 = diff([31440.7 31440.7]);
        dist23 = sqrt(sum([dX23 dY23 dZ23].^2));
        % Distance between beads 1 and 3
        dX13 = diff([309 221]);
        dY13 = diff([80 -22]);
        dZ13 = diff([31440.7 31440.7]);
        dist13 = sqrt(sum([dX13 dY13 dZ13].^2));
        % ========BEAD 1
        % Double bead; 2/5 of BM
        % CF = 20.08
        % Indices to groups of recs
        B1iBB = {[1 43:45] [2 46:48] [3 49:51] [4 52:53] [5 54:55] 73 75 77}; % 0.1-30 kHz
        B1i11_22 = {56:58 59:61 62:64 65:67 68:70 71:72 74 76 78}; % 11k-22k5
        % apples
        AP_1w = apple(DS(B1iBB),2,0,ME)
        AP_1n = apple(DS(B1i11_22),2,0,ME)
        % ========BEAD 2
        % Solitary bead; 2/5 of BM
        % CF = 18.65
        % Indices to groups of recs
        B2iBB = {6:9 10:13 14:16 17:20 21:24 79 81 83}; % 0.1-30 kHz
        B2i11_22 = {28:30 25:27 31:33 34:36 37:39 40:42 80 82 84}; % 11k-22k5
        % apples
        AP_2w = apple(DS(B2iBB),2,0,ME)
        AP_2n = apple(DS(B2i11_22),2,0,ME)
        % ========BEAD 3
        % Solitary bead; 1/4 of BM
        % bead landed during measurements of previous beads
        % CF = 17.06
        % Indices to groups of recs
        B3iBB = {86 89:90 93:94 97:98 101:102 105 107 109}; % 0.1-30 kHz
        B3i11_22 = {87:88 91:92 95:96 99:100 103:104 106 108 110}; % 11k-22k5
        % apples
        AP_3w = apple(DS(B3iBB),2,0,ME)
        AP_3n = apple(DS(B3i11_22),2,0,ME)
        % ==============TRANSFER FUNCTIONS==========
        switch upper(ExpName),
            case 'RG12438_12' % bead 1 -> 2
                dX = dX12;
                dY = dY12;
                dZ = dZ12;
                BeadDist = dist12;
                AP_Aw = AP_1w;
                AP_An = AP_1n;
                AP_Bw = AP_2w;
                AP_Bn = AP_2n;
                Lw = CVlocalTRF(DS(B1iBB),DS(B2iBB),2,BeadDist,noPlot) % 0.1-30 kHz
                Ln =CVlocalTRF(DS(B1i11_22),DS(B2i11_22),2,BeadDist,noPlot) % 11k-22k5
            case 'RG12438_23' % bead 2 -> 3
                dX = dX23;
                dY = dY23;
                dZ = dZ23;
                BeadDist = dist23;
                AP_Aw = AP_2w;
                AP_An = AP_2n;
                AP_Bw = AP_3w;
                AP_Bn = AP_3n;
                Lw = CVlocalTRF(DS(B2iBB),DS(B3iBB),2,BeadDist,noPlot) % 0.1-30 kHz
                Ln = CVlocalTRF(DS(B2i11_22),DS(B3i11_22),2,BeadDist,noPlot) % 11k-22k5
            case 'RG12438_13' % bead 1 -> 3
                dX = dX13;
                dY = dY13;
                dZ = dZ13;
                BeadDist = dist13;
                AP_Aw = AP_1w;
                AP_An = AP_1n;
                AP_Bw = AP_3w;
                AP_Bn = AP_3n;
                Lw = CVlocalTRF(DS(B1iBB),DS(B3iBB),2,BeadDist,noPlot) % 0.1-30 kHz
                Ln = CVlocalTRF(DS(B1i11_22),DS(B3i11_22),2,BeadDist,noPlot) % 11k-22k5
        end
    otherwise,
        error(['Unknown experiment ' ExpName '.']);
end % switch upper(EXP)
if ~exist('LnNames', 'var'), LnNames = {'Ln'}; end
local_strip_APs;
AP_Adum = []; AP_Bdum = []; % avoid crash in next line
save(Sfile, 'ExpName', 'L*', 'dX', 'dY', 'dZ', 'BeadDist', 'AP_A*', 'AP_B*',  'ME');


function AP = local_strip_anAP(AP);
[AP.Reference] = deal([]);

function local_strip_APs;
qq = evalin('caller', 'who(''AP_*'')')
for ii=1:numel(qq),
    cmd = [qq{ii} ' = local_strip_anAP(' qq{ii}  ');'];
    evalin('caller', cmd);
end


function Y = local_get(ExpName);
Sfile = fullfile(local_Sdir, [ExpName '_localTRF.mat']);
Y = load(Sfile);
% add integer # cycles to some phase curves
switch upper(ExpName),
    case 'RG12401',
        Y.Lw.A2B(3).Phase = Y.Lw.A2B(3).Phase+1 ;
        Y.Lw.A2B(9).Phase = Y.Lw.A2B(9).Phase+1 ;
        Y.L3.A2B(1).Phase = Y.L3.A2B(1).Phase+1 ;
        Y.L3.A2B(2).Phase = Y.L3.A2B(2).Phase + 1;
        Y.L3.A2B(3).Phase = Y.L3.A2B(3).Phase + 1;
        Y.L15.A2B(1).Phase = Y.L15.A2B(1).Phase - 1 ;
        Y.L15.A2B(2).Phase = Y.L15.A2B(2).Phase - 1;
        Y.L15.A2B(4).Phase = Y.L15.A2B(4).Phase + 1;
        Y.L15.A2B(5).Phase = Y.L15.A2B(5).Phase + 1;
    case 'RG12405',
        for jj=6:9, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=1:9, Y.L1.A2B(jj).Phase = Y.L1.A2B(jj).Phase + 1; end
        for jj=1:6, Y.L15.A2B(jj).Phase = Y.L15.A2B(jj).Phase - 1; end
        for jj=8:9, Y.L15.A2B(jj).Phase = Y.L15.A2B(jj).Phase +  1; end
    case 'RG12407',
        for jj=3:9, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=1:9, Y.L1.A2B(jj).Phase = Y.L1.A2B(jj).Phase + 1; end
        for jj=1:9, Y.L3.A2B(jj).Phase = Y.L3.A2B(jj).Phase + 1; end
        for jj=1:4, Y.L6.A2B(jj).Phase = Y.L6.A2B(jj).Phase + 1; end
        for jj=9, Y.L15.A2B(jj).Phase = Y.L15.A2B(jj).Phase + 1; end
    case 'RG12417',
        for jj=4:9, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=5, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=4:9, Y.L1.A2B(jj).Phase = Y.L1.A2B(jj).Phase + 1; end
        %for jj=4:9, Y.L3.A2B(jj).Phase = Y.L1.A2B(jj).Phase + 1; end
    case 'RG12422',
        for jj=4:9, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        Y.Ln = {Y.L1 Y.L4 Y.L10 Y.L16 Y.L22};
        for jj=setdiff(1:9,6), Y.L1.A2B(jj).Phase = Y.L1.A2B(jj).Phase + 1; end
        for jj=[4 6 7 8 9], Y.L4.A2B(jj).Phase = Y.L4.A2B(jj).Phase + 1; end
        for jj=[8 9], Y.L16.A2B(jj).Phase = Y.L16.A2B(jj).Phase + 1; end
        for jj=1:3, Y.L22.A2B(jj).Phase = Y.L22.A2B(jj).Phase - 1; end
        for jj=4:8, Y.L22.A2B(jj).Phase = Y.L22.A2B(jj).Phase + 1; end
    case 'RG12423',
        for jj=1:7, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=[3 4 6], Y.L1.A2B(jj).Phase = Y.L1.A2B(jj).Phase + 1; end
        for jj=3, Y.L4.A2B(jj).Phase = Y.L4.A2B(jj).Phase + 1; end
        for jj=[4 5 6], Y.L16.A2B(jj).Phase = Y.L16.A2B(jj).Phase + 1; end
        for jj=[5 6 7], Y.L22.A2B(jj).Phase = Y.L22.A2B(jj).Phase + 1; end
    case 'RG12424',
        for jj=3:9, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=1:9, Y.L1.A2B(jj).Phase = Y.L1.A2B(jj).Phase + 1; end
        for jj=[4 6:9], Y.L4.A2B(jj).Phase = Y.L4.A2B(jj).Phase + 1; end
        for jj=5:8, Y.L22.A2B(jj).Phase = Y.L22.A2B(jj).Phase + 1; end
        for jj=1:4, Y.LHFP.A2B(jj).Phase = Y.LHFP.A2B(jj).Phase + 1; end
    case 'RG12432',
        for jj=[3 5], Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
    case 'RG12433',
        for jj=6:8, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=7:8, Y.Ln.A2B(jj).Phase = Y.Ln.A2B(jj).Phase + 1; end
    case 'RG12434',
        for jj=2:8, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        %for jj=7:8, Y.Ln.A2B(jj).Phase = Y.Ln.A2B(jj).Phase + 1; end
    case 'RG12436',
        for jj=[2 4:8], Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        %for jj=7:8, Y.Ln.A2B(jj).Phase = Y.Ln.A2B(jj).Phase + 1; end
    case 'RG12438_12',
        for jj=1:8, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=9, Y.Ln.A2B(jj).Phase = Y.Ln.A2B(jj).Phase + 1; end
    case 'RG12438_23',
        for jj=1:8, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=7:8, Y.Ln.A2B(jj).Phase = Y.Ln.A2B(jj).Phase + 1; end
    case 'RG12438_13',
        for jj=1:8, Y.Lw.A2B(jj).Phase = Y.Lw.A2B(jj).Phase + 1; end
        for jj=8, Y.Ln.A2B(jj).Phase = Y.Ln.A2B(jj).Phase + 1; end
    otherwise,
        warning('No phase correction was prepared for this experiment.');
end % switch upper(ExpName)
Lnames = [{'Lw'} Y.LnNames];
for jj=1:numel(Lnames),
    Y.(Lnames{jj}) = local_c(Y.(Lnames{jj}));
    Y.(Lnames{jj}) = local_U(3,1,Y.(Lnames{jj}));
    Y.(Lnames{jj}) = local_addMid(Y.(Lnames{jj}));
end

function varargout = local_c(varargin);
% recompute phase velocity, now using correct phases
% also supply wavenumber k (cycle/m) and wavelength (um)
for iarg=1:nargin,
    L = varargin{iarg};
    for ii=1:numel(L.A2B),
        L.A2B(ii).c = -1e-6*L.A2B(ii).DX.*L.A2B(ii).Fprim./L.A2B(ii).Phase; % m/s
        L.A2B(ii).k = L.A2B(ii).Fprim./L.A2B(ii).c; % cycle/m
        L.A2B(ii).lambda = 1e6./L.A2B(ii).k;
    end
    L.A2B = structmovefield(L.A2B, 'k', 'cD');
    L.A2B = structmovefield(L.A2B, 'lambda', 'cD');
    varargout{iarg} = L;
end

function varargout = local_U(Nfit, fitOrd, varargin);
for iarg=1:nargin-2,
    L = varargin{iarg};
    for ii=1:numel(L.A2B),
        ab = L.A2B(ii);
        freq = ab.Fprim/1e3;
        tau = ab.DX./ab.c;
        tauge = local_slope(freq, freq.*tau, Nfit, fitOrd);
        L.A2B(ii).Uo = L.A2B(ii).U; 
        L.A2B(ii).tau = tau; % phase delay in us
        L.A2B(ii).tauG = tauge; % group delay in us
        L.A2B(ii).U = L.A2B(ii).DX./tauge; % group velocity in m/s
        L.A2B(ii).U_over_c = L.A2B(ii).U./L.A2B(ii).c;
    end
    L.A2B = structmovefield(L.A2B, 'U', 'cD');
    L.A2B = structmovefield(L.A2B, 'tau', 'cD');
    L.A2B = structmovefield(L.A2B, 'tauG', 'cD');
    varargout{iarg} = L;
end

function S = local_slope(x,y,Np, fitOrd);
Np = 1+2*floor(Np/2);
N2 = floor(Np/2);
S = nan(size(x));
N = numel(x);
for ii=1:N,
    ir = ii+(-N2:N2);
    if min(ir)<1, ir = 1+ir-min(ir); end
    if max(ir)>N, ir = N+ir-max(ir); end
    xf = x(ir); yf = y(ir);
    [yf, dum, xf] = denan(yf, xf);
    if numel(xf)<3,
        S(ii) = nan;
    else,
        Xm = mean(xf);
        P = polyfit(xf-Xm, yf, fitOrd);
        Pd = PolyDiff(P);
        Xi = x(ii)-Xm;
        S(ii) = polyval(Pd, x(ii)-Xm);
        Si = S(ii);
    end
end

function varargout = local_addMid(varargin);
for ia=1:nargin,
    L = varargin{ia};
    for ii=1:numel(L.AP_A),
        apa = L.AP_A(ii);
        apb = L.AP_B(ii);
        apm = apa;
        apm.Gain = 0.5*(apa.Gain+apb.Gain);
        apm.Phase = 0.5*(apa.Phase+apb.Phase);
        apm.Alpha = max(apa.Alpha, apb.Alpha);
        L.AP_mid(ii) = apm;
    end
    varargout{ia} = L;
end

function [Ln, AP_An, AP_Bn] = local_Ln(Y);
% narrowband TRFs as cell array
N = numel(Y.LnNames);
[AP_An, AP_Bn] = deal({});
for ii=1:N,
    lname = Y.LnNames{ii};
    Ln{ii} = Y.(lname);
    aname = strrep(lname, 'L', 'AP_A');
    if isfield(Y, aname), 
        AP_An = [AP_An Y.(aname)];
    end
    bname = strrep(lname, 'L', 'AP_B');
    if isfield(Y, bname), 
        AP_Bn = [AP_Bn Y.(bname)];
    end
end


function local_TRFplot(L, ha);
% plot transfer functions
L = cellify(L);
for ii=1:numel(L),
    CVlocalTRF(L{ii}, 'plotns', ha);
end
axes(ha(1)); 
legend off;
xlog125([1 40]);
ylabel('Gain (dB)');
axes(ha(2));
ylim([-3.5 0.5]);
xlabel('Frequency (kHz)');
ylabel('Phase (cycle)');

function local_dispersion_plot(L, ha);
% plot frequency vs wavenumber
axes(ha);
L = cellify(L);
kmax = 0;
for iL=1:numel(L),
    Li = L{iL};
    Nw = numel(Li.A2B);
    for ii=Nw:-1:1,
        ab = Li.A2B(ii);
        xplot(ab.k/1e3, ab.Fprim/1e3+pmask(ab.Alpha<=0.001), jetcol(ii/(Nw+1)), 'linewidth', 2);
        kmax = max(kmax, max(ab.k/1e3+pmask(ab.Alpha<=0.001)));
    end
end
xlim(kmax*[-0.1 1.1]);
ylim([0 1.1*max(ab.Fprim)/1e3]);
grid on;
xlabel('Wavenumber (cycles/mm)');
ylabel('Frequency (kHz)');

function local_GainPhase_plot(L, ha);
% plot gain vs phase
axes(ha);
L = cellify(L);
Cph = A2dB(exp(2*pi)); % constant relating 1 cycle to =54.4 dB
phmin = inf;
for iL=1:numel(L),
    Li = L{iL};
    Nw = numel(Li.A2B);
    PH = []; GA = [];
    for ii=Nw:-1:1,
        ab = Li.A2B(ii);
        ph = -Cph*ab.Phase;
        ga = ab.Gain+pmask(ab.Alpha<=0.001);
        xplot(-Cph*ab.Phase, ga, jetcol(ii/(Nw+1)), 'linewidth', 2);
        phmin = min(phmin, min(ab.Phase+pmask(ab.Alpha<=0.001)));
        PH = [PH; ph];
        GA = [GA; ga];
    end
    xplot(PH, GA, 'k')
    set(gca, 'xtick', -1000:20:1000, 'ytick', -1000:20:1000);
    xlim(Cph*[-0.25 -1.1*phmin]);
    axis equal;
    grid on
    xlabel('Phase (cycle/54.4)');
    ylabel('Gain (dB)');
end

function local_Velo_plot(Velotype, L, ha);
% plot velocity vs freq
if isequal('phase', Velotype), vfield = 'c';
elseif isequal('group', Velotype), vfield = 'U';
else, error(VeloType);
end
axes(ha);
L = cellify(L);
for iL=1:numel(L),
    Li = L{iL};
    Nw = numel(Li.A2B);
    for ii=Nw:-1:1,
        ab = Li.A2B(ii);
        xplot(ab.Fprim/1e3, ab.(vfield)+pmask(ab.Alpha<=0.001), jetcol(ii/(Nw+1)), 'linewidth', 2);
    end
end
xlog125([1 40]);
xlabel('Frequency (kHz)');
if isequal('phase', Velotype), 
    ylabel('Phase velocity (m/s)');
else
    ylabel('Group velocity (m/s)');
end
xlog125([1 40]);
ylog125([0.5 200]);
grid on

function local_Group_over_Phase_plot(L, ha);
% plot group velocity vs phase velocity
axes(ha);
L = cellify(L);
for iL=1:numel(L),
    Li = L{iL};
    Nw = numel(Li.A2B);
    for ii=Nw:-1:1,
        ab = Li.A2B(ii);
        xplot(ab.Fprim/1e3, ab.U./ab.c+pmask(ab.Alpha<=0.001), jetcol(ii/(Nw+1)) );
        %freq = ab.Fprim/1e3;
        %ffreq = logispace(min(freq), max(freq), 1e3);
        %k = interp1(freq, ab.k, ffreq);
        %Alpha = interp1(freq, ab.Alpha, ffreq);
        %U_c = diff(log(ffreq))./diff(log(k));
        %xplot(ffreq(2:end), U_c+pmask(Alpha(2:end)<=0.001), jetcol(ii/(Nw+1)) );
    end
end
xlabel('Frequency (kHz)');
ylabel('U/c');
xlog125([1 40]);
ylog125([0.1 2]);
grid on

function local_Group_vs_Phase_plot(L, ha);
% plot group velocity vs phase velocity
axes(ha);
L = cellify(L);
for iL=1:numel(L),
    Li = L{iL};
    Nw = numel(Li.A2B);
    for ii=Nw:-1:1,
        ab = Li.A2B(ii);
        xplot(ab.U, ab.c+pmask(ab.Alpha<=0.001), jetcol(ii/(Nw+1)) );
    end
end
xlabel('Phase velocity (m/s)');
ylabel('Group velocity (m/s)');
xlog125([0.5 200]);
ylog125([0.5 200]);
grid on

function local_lambda_plot(L, ha);
% plot group velocity vs phase velocity
axes(ha);
L = cellify(L);
for iL=1:numel(L),
    Li = L{iL};
    Nw = numel(Li.A2B);
    for ii=Nw:-1:1,
        ab = Li.A2B(ii);
        xplot(ab.Fprim/1e3, ab.lambda/1e3+pmask(ab.Alpha<=0.001), jetcol(ii/(Nw+1)), 'linewidth',2);
    end
end
xlabel('Frequency (kHz)');
ylabel('Wavelength (mm)');
xlog125([0.5 40]);
ylog125([0.05 50]);
grid on

function local_Group_vs_Lambda_plot(L, ha);
% plot group velocity vs lambda
axes(ha);
L = cellify(L);
for iL=1:numel(L),
    Li = L{iL};
    Nw = numel(Li.A2B);
    for ii=Nw:-1:1,
        ab = Li.A2B(ii);
        xplot(ab.lambda/1e3, ab.U+pmask(ab.Alpha<=0.001), jetcol(ii/(Nw+1)), 'linewidth', 2);
    end
end
xlabel('Wavelength (mm)');
ylabel('Group velocity (m/s)');
xlog125([0.05 20]);
ylog125([0.5 200]);
grid on


function LH = local_plot_single(AP, ha, LW);
% plot stapes-BM gain & phase of single bead
AP = cellify(AP);
for iL=1:numel(AP),
    APi = AP{iL};
    LH = apple(dataset(), APi, 'plot', ha);
    delete([LH.Gain_dot LH.Gain_dot]);
    set([LH.Gain_lin LH.Phase_lin], 'markersize', 3, 'linewidth', LW);
end
axes(ha(1));
legend off;
xlog125([0.4 40]);
axes(ha(2));
legend off;
xlog125([0.4 40]);
linkaxes(ha, 'x');

function local_plot1(L, AP_A, AP_B, T);
% ===plot first set of analyses: local transfers, resp re stapes
for ii=1:4, ha(ii) = subplot(2,2,ii); end
local_TRFplot(L, ha([1 3]));
local_plot_single(AP_A, ha([2 4]), 2);
local_plot_single(AP_B, ha([2 4]), 1);
title(ha(1), T, 'fontsize', 14);

function local_plot2(L, T);
% ===plot 2nd set of analyses: dispersion diagram, Gain vs Phase, U & c
for ii=1:8, ha(ii) = subplot(2,4,ii); end
% dispersion diagram
local_dispersion_plot(L, ha(1));
axes(ha(1));
title(T, 'fontsize', 14);
% Gain versus amplitude
local_GainPhase_plot(L, ha(5));
% Phase and group velocity
local_Velo_plot('phase', L, ha(2));
local_Velo_plot('group', L, ha(6));
% group over phase velocity
local_Group_over_Phase_plot(L, ha(3));
% group vs phase velocity
local_Group_vs_Phase_plot(L, ha(7));
% wavelength vs freq
local_lambda_plot(L, ha(4));
% group vs k
local_Group_vs_Lambda_plot(L, ha(8));

function local_show(ExpName);
Y = local_get(ExpName);
TT = sprintf('%s;   Bead Distance %d \\mum', Y.ExpName, round(Y.BeadDist));
Sdir = fullfile(processed_datadir, 'CV', 'twobeads', 'figs');
% ==Wideband
if ~isfield(Y, 'AP_Aw'),
    Y.AP_Aw = {};
    Y.AP_Bw = {};
end
fh1 = figure;
set(fh1,'units', 'normalized', 'position', [0.0141 0.16 0.489 0.745])
local_plot1(Y.Lw, Y.AP_Aw, Y.AP_Bw, TT);
saveas(fh1, fullfile(Sdir, [ExpName '_WB_1.fig']));
fh2 = figure;
set(fh2,'units', 'normalized', 'position', [0.01 0.15 0.964 0.745]);
local_plot2(Y.Lw, TT); 
saveas(fh2, fullfile(Sdir, [ExpName '_WB_2.fig']));
% ==Narrowband
[Ln, AP_An, AP_Bn] = local_Ln(Y);
fh3 = figure;
set(fh3,'units', 'normalized', 'position', [0.034 0.18 0.489 0.745])
local_plot1(Ln, AP_An, AP_Bn, TT);
saveas(fh3, fullfile(Sdir, [ExpName '_NB_1.fig']));
fh4 = figure;
set(fh4,'units', 'normalized', 'position', [0.03 0.17 0.964 0.745]);
local_plot2(Ln, TT); 
saveas(fh4, fullfile(Sdir, [ExpName '_NB_2.fig']));

function local_cfmatch(ExpName, dA, FreqRat);
set(figure,'units', 'normalized', 'position', [0.23 0.0611 0.518 0.853]);
Y = local_get(ExpName);
for ii=1:2, ha(ii) = subplot(2,1,ii); end;
for ii=1:numel(Y.AP_Aw),
    Y.AP_Aw(ii).Fprim = Y.AP_Aw(ii).Fprim/sqrt(FreqRat);
    Y.AP_Aw(ii).Gain = Y.AP_Aw(ii).Gain-dA/2;
end
for ii=1:numel(Y.AP_Bw),
    Y.AP_Bw(ii).Fprim = Y.AP_Bw(ii).Fprim*sqrt(FreqRat);
    Y.AP_Bw(ii).Gain = Y.AP_Bw(ii).Gain+dA/2;
end
LHA = local_plot_single(Y.AP_Aw, ha, 2);
LHB = local_plot_single(Y.AP_Bw, ha, 1);
SPL = [Y.AP_Aw.baseSPL Y.AP_Bw.baseSPL];
LH = [LHA.Gain_lin, LHB.Gain_lin];
[dum, dum, iun] = unique(SPL);
for ii=1:numel(LH),
    set(LH(ii), lico(iun(ii)));
end
set(ha, 'xlim', [5 28], 'xtick', 1:30, 'xticklabel', 1:30);
grid(ha(1), 'on'); grid(ha(2), 'on'); 









