function [D, fac, CMgain] = correctCMamplitude(D);
%correctCMamplitude -correct amplitude in case of CM recordings
%
%D = correctCMamplitude(D) "converts" the amplitude data in struct D (from ZWOAEimport)
% from "dB SPL" to dB re. 1V. This is for CM recordings. The CM
% amplification factor is also taken into account.
%
% [D, fac] = correctCMamplitude(D) also returns the factor fac (in dB) that
% has to be added to the amplitude data to go from dB SPL to dBV.
% This factor includes both dB SPL --> dBV and amplification
%
% [D, fac, ampl] = correctCMamplitude(D) also returns the conversion factor
% ampl (which is already included in fac!). This is the correction for only the amplification (in dB)
%
% NOTE: this function is already invoked within ZWOAEimport so that CM
% recordings are "corrected" when returned by ZWOAEimport.
%
%

fac = 0;
ampl = 0;

%do only when recording is 'CM'
if strcmp(D.RecType,'CM') == 1,
    fac = D.micgain+D.dBV-D.dBSPL;     %from dB SPL to dB re 1V
        
    %also correct for CM amplification
    [FN PN] = ZWOAEfilename(D.ExpID, D.recID,[],'raw');
    load(fullfile(PN, FN),'-mat'); %all data/info loaded into struct ST
    tmp = ST.Setupinfo.CMgain; %get the CM amplification settings
    CMgain = -sum(A2dB(tmp(~isinf(tmp)))); % .CMgain may be an array with gain settings for consecutive amplifiers
    fac = fac+CMgain;
    
end

D.MG = D.MG+fac; %apply the correction factor to amplitude spectrum

    