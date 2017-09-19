function Y = CopySelectedTKdata;

DestDir = 'D:\USR\Marcel\sortTKdata\tempABFstore';
Sall = load('D:\USR\Marcel\sortTKdata\all.SGSR_ABFs', '-mat'); Sall=Sall.S;

local_get_oneSGSR('initialize', DestDir); % initialize local fcn
%-----RG09110 cell 3 -----05-May-2009 16:06:02-
ExpID = 'RG09110';
%[a b FullABFdir c] = SGSR2ABF(Sall, ExpID, '
% '); 
%
local_get_oneSGSR(Sall, ExpID, '3-1-FS'); 
local_get_oneSGSR(Sall, ExpID, '3-2-FS'); 
local_get_oneSGSR(Sall, ExpID, '3-3-FS'); 
local_get_oneSGSR(Sall, ExpID, '3-5-FS'); 
local_get_oneSGSR(Sall, ExpID, '3-6-FS'); 
local_get_oneSGSR(Sall, ExpID, '3-7-FS'); 
local_get_oneSGSR(Sall, ExpID, '3-8-FS'); 
local_get_oneSGSR(Sall, ExpID, '3-9-FS'); 
local_get_oneSGSR(Sall, ExpID, '3-10-FS'); 
local_get_oneSGSR(Sall, ExpID, '3-15-SPL'); 
local_get_oneSGSR(Sall, ExpID, '3-16-SPL'); 
local_get_oneSGSR(Sall, ExpID, '3-33-SPL'); 
local_get_oneSGSR(Sall, ExpID, '3-34-SPL'); 

% ======== COMPLEX A ======================
%    - CF ~ 2 kHz
%    - tau = 3.5 ms (BN 14)
%    - nonmonotonic rate fnc
%    - *** prepotential
% 240409pen10u13_RG09106\recording104__14_22_SPL
% 240409pen10u13_RG09106\recording100__14_18_FS

% ------RG09106 cell 14 -----24-Apr-2009 23:39:13-----
ExpID = 'RG09106';
local_get_oneSGSR(Sall, ExpID, '14-1-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-2-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-3-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-4-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-5-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-6-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-7-SPL'); 
local_get_oneSGSR(Sall, ExpID, '14-10-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-11-PS'); 
local_get_oneSGSR(Sall, ExpID, '14-13-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-14-SPL'); 
local_get_oneSGSR(Sall, ExpID, '14-15-SPL'); 
local_get_oneSGSR(Sall, ExpID, '14-16-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-17-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-18-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-19-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-20-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-21-FS'); 
local_get_oneSGSR(Sall, ExpID, '14-22-SPL'); 


% ======== Complex B ==========
%    - CF ~ 1100 Hz (BN 1,2)
%    - highly nonmonotonic rate fcn
%    - Tau = 4.5 ms (BN 1,2)
% 
% 020609pen2u1_wtg002\recording2__1_2_FS
% 020609pen2u1_wtg002\recording19__1_28_SPL
% 
% ------WTG002 cell 1 -----02-Jun-2009 14:22:31-----
ExpID = 'WTG002';
local_get_oneSGSR(Sall, ExpID, '1-1-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-2-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-3-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-4-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-5-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-6-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-7-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-11-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-13-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-14-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-15-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-16-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-17-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-20-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-24-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-25-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-26-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-27-FS'); 
local_get_oneSGSR(Sall, ExpID, '1-28-SPL'); 
local_get_oneSGSR(Sall, ExpID, '1-29-SPL'); 
local_get_oneSGSR(Sall, ExpID, '1-30-SPL'); 


% ======== Low CF ==========
% 
%    - CF ~ 650 Hz (BN 3)
%    - monotonic rate fcn
%    - tau = 5.3 ms (BN 3)
% 
% 020609pen2u2_wtg002\recording21__2_1_FS
% 020609pen2u2_wtg002\recording35__2_19_SPL
% 
% ------WTG002 cell 2 -----02-Jun-2009 15:26:30
local_get_oneSGSR(Sall, ExpID, '2-1-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-2-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-3-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-4-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-5-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-6-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-7-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-12-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-13-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-14-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-14-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-15-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-16-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-17-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-18-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-19-SPL'); 
local_get_oneSGSR(Sall, ExpID, '2-20-SPL'); 
local_get_oneSGSR(Sall, ExpID, '2-22-SPL'); 
local_get_oneSGSR(Sall, ExpID, '2-25-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-26-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-27-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-28-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-29-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-30-FS'); 
local_get_oneSGSR(Sall, ExpID, '2-31-FS'); 


% ======== High_Sync ==========
% 
%    - CF ~ 300 Hz
%    - tau = 7 ms (BN 8,9,10)
%    - monotonic rate fcn
%    - high sync
%    - no prepotentials
%    - no clear EPSPs
% 
% 020609pen4u4_wtg002_highsync\recording64__4_8_FS
% 020609pen4u4_wtg002_highsync\recording65__4_11_FS
% 020609pen4u4_wtg002_highsync\recording66__4_15_SPL

% ------WTG002 cell 4 -----02-Jun-2009 19:31:41-----el(1):NaN um---
local_get_oneSGSR(Sall, ExpID, '4-1-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-2-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-3-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-4-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-5-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-8-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-11-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-15-SPL'); 
local_get_oneSGSR(Sall, ExpID, '4-18-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-19-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-20-FS'); 
local_get_oneSGSR(Sall, ExpID, '4-21-FS'); 

Y = local_get_oneSGSR('getstock');
%=========================
function Y=local_get_oneSGSR(Sall, ExpID, recID);
persistent S_stock DestDir
if isequal('initialize', Sall), % local_get_oneSGSR(Sall, ExpID, 'initialize', DestDir);
    S_stock = [];
    DestDir = ExpID;
    return;
elseif isequal('getstock', Sall), % 
    Y = S_stock;
    return;
end
S = SGSR_ABF_filter(Sall, ExpID, recID);
disp('============')
[ABFdir, ABFfiles, FullABFdir] = ABF_fileInfo(S)
%copyfile(fullfile(FullABFdir, '*.ABF'), DestDir);
S_stock = [S_stock, S];





