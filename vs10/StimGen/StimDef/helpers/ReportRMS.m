function RMSstr = ReportRMS(figh, P)
% ReportFsam - report rms value xRMS in RMSDisp panel 
%     ReportRMS(figh, xRMS) reports the rms value xRMS (in
%     mV) to the RSMDisp messenger of the stimulus GUI having handle figh.
%
%     ReportFsam(figh, nan) resets the string of the FsamDisp messenger 
%     to its TestLine value (see Messenger).
%
%   See StimGUI, PlayTime, Messenger.

% compute and format
    RMSstr = {'', ''};

    if isfield(P,'RMSav')
        xRMS = P.RMSav*1e3;
        DAC = P.DAC;
        nChan = numel(xRMS);
        str = '%s rms = %0.1f mV';
        
        if nChan < 2
            DAC = cellstr(DAC);
        else
            DAC = {'R', 'L'}; % to be Confirmed
        end
    for ch=1:nChan
        RMSstr{ch} = sprintf(str, DAC{ch}(:), xRMS(ch));
    end
    end

% report
for i=1:2
M = find(messenger(), figh, ['RMSDisp' num2sstr(i)]);
report(M,RMSstr{i});
end



