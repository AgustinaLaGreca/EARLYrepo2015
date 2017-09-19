% -------------------------------------------------------------------------
% Getting the latest EARLY from the server 
% -------------------------------------------------------------------------
% There should be 3 folders in the
%   C: drive:
% 	C:\EARLY                    (general installation folder)
% 	C:\Early_StimDefLeuven      (stim definitions) 
%    C:\Early_Systemdata        (setup specific files)
% 	
% 
% To reset to the latest version on the server Remove all files/folders
% from C:\EARLY except for a folder called “.svn”. The later is normally 
% hidden so don’t worry if you can’t find it.
% 
% Update the C:\EARLY folder by right clicking on it, select “Update” from
% the “tortoiseSVN” menu. If the system asks for a password use:
% user=Experiment/password=1monkey.
% 
% Follow the same procedure for the Early_StimDefLeuven folder
% 
% If the C:\Early_Systemdata is missing, just copy it back from
% EARLY\CpuSpecificFiles\Early_Systemdata\
% 
% -------------------------------------------------------------------------	
% Reverting to a previous EARLY version 
% -------------------------------------------------------------------------	
% If the problem persists after you followed the procedure  above, 
% there might be a problem with the actual
% version on the server. A useful tool to determine to which version you
% should revert is the logViewer from the tortoiseSVN menu. It  shows the
% comments entered by the developer at each version change.
% 
% To revert to a previous version, using the following procedure: Right
% click the EARLY folder and select “Log Viewer” from the  tortoiseSVN
% menu. Select the version you want to revert to and right click on the
% line within the Log Viewer window Select revert to this version from the
% menu
% 
