function out = ZWOAEdatadir(PN)
% ZWOAEdatadir -get/set directory for saving ZWOAE data
%
%   out = ZWOAEdatadir returns the data save directory for ZWOAE as is
%   currently specified within the setupfile ('DataRootDirs'). If not specified, an error is
%   cast. out is the returned directory
%
%   out = ZWOAEdatadir(PN) sets the data save directory for ZWOAE to PN. If the specified
%   directory does not exist, it is created. out is this new directory
%
%   Syntax examples:
%       DD = ZWOAEdatadir;
%       DD = ZWOAEdatadir('d:\data\ZWOAE');
%
%   See also ToSetupFile, FromSetupFile, SetupList


SUF = 'DataRootDirs'; % local setup file for datadir listing
w = SetupList; %get all setup files 

if nargin == 0,
    %no input, return currently saved datadir for ZWOAE

    % check whether SUF exists within all setup files. If not, cast error
    [w, mess] = keywordMatch(SUF, w);
    if ~isempty(mess) == 1,
        error('No setupfile for DataRootDirs');
    end
    out = FromSetupFile(SUF, 'ZWOAE');
elseif nargin == 1 && ischar(PN) == 1,
    % directory must exist
    if ~exist(PN,'dir'),
        error(['Non existent directory ''' PN '''.']);
    end
    ToSetupFile(SUF, '-propval', 'ZWOAE', PN);
    out = FromSetupFile(SUF, 'ZWOAE');
    
else
    error('Input must be a string specifying absolute path');
end

