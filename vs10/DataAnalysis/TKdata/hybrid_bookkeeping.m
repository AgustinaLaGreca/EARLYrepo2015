function recording_infos = hybrid_bookkeeping(DatasetDir)
% hybrid_bookkeeping - link data from SGSR and Clampex
%   S = hybrid_bookkeeping(D) looks for .abf files in directory D and links
%   them to SGSR datasets by using Thomas' "log files". It also retrieves
%   additional info on the state of the Multiclamp amplifier from Thomas'
%   "state files". The following conventions are needed for a correct
%   retrieval of this info:
%     1. The grandparent dir of D is the "Experiment directory" whose 
%        sub-sub directories toghether contain all the ABF files of the 
%        experiment. The Experiment Directory is called ExpDir below.
%     2. ExpDir contains a subdir <ExpDir>\logs containing all the log
%        files of the experiment.
%     3. ExpDir contains a subdir <ExpDir>\statefiles containing all the 
%        state files of the experiment.
%      
%  Hybrid_bookkeeping returns a struct S containing the info on the 
%  SGSR dataset corresponding to the ABF files in directory D, and the info
%  from the most recent state file.
%
%  Example:
%
%    hybrid_bookkeeping('D:\USR\Jeannette\data\RG09163\unit7\record1')
%    
% ans = 
%        Date_of_recording: '07/10/09'
%      Timestamp_first_abf: '22:34:46'
%       Timestamp_last_abf: '22:35:18'
%                   d_____: '__________'
%                SGSRExpID: 'RG09163'
%           SGSRidentifier: '6-1-FS'
%       Stim_Property_Unit: 'Hz'
%     Stim_Property_Values: [25x1 double]
%              Rec_Quality: 'Membrane Potential'
%            Rec_Base_unit: '10 V/V'
%                 Rec_Gain: 10
%
%  See also compile_all_SGSR_ABFs, TKpool.

if ~exist(DatasetDir, 'dir'), %error checking
    error(['Directory ''',DatasetDir, ''' does not exist! Try "dir" or "win-key + e"...']);
end

ExpDir = fileparts(fileparts(DatasetDir)); % grandparent dir = experiment dir
abffiles = dir(fullfile(DatasetDir, '*.abf'));
%[log_path, state_path] = parse_infopaths(PATH);
%state_path =
log_path = fullfile(ExpDir, 'logs');
logfiles = dir(fullfile(log_path, '*.log'));
if isempty(logfiles),
    error(['No log files found in ''' log_path '''.']);
end
state_path = fullfile(ExpDir, 'statefiles');
statefiles = dir(fullfile(state_path, '*.state'));
if isempty(statefiles),
    error(['No state files found in ''' state_path '''.']);
end

start_rec = min([abffiles.datenum]); % datenum of oldest abf, i.e., the first recorded
stop_rec = max([abffiles.datenum]); % datenum of newest abf, i.e., the last recorded
%rec_time = stop_rec - start_rec;
%rec_time = datevec(rec_time);
%disp(['These recordings took ' num2str(rec_time(5)) ' min ' num2str(rec_time(6)) 'sec']);
if isempty(start_rec),
    recording_infos = [];
    return;
end
%find & read the youngest logfile that is older than our recording
q_older = [logfiles.datenum] < start_rec; % logical indexing logfiles, indicating older than start_rec
oldlogfiles = logfiles(q_older); % selection of older logs
[dum, imax] = max([oldlogfiles.datenum]); 
logname = [oldlogfiles(imax).name]; % name of targeted log file
log_head = logname; % store for parsing the name
%log_path, logname
log_content = mytextread(fullfile(log_path, logname));

%find the youngest statefile that is older than our recording
q_older = [statefiles.datenum] < start_rec; % logical indexing statefiles, indicating older than start_rec
oldstatefiles = statefiles(q_older); % selection of older states
[dum, imax] = max([oldstatefiles.datenum]); 
statename = [oldstatefiles(imax).name]; % name of targeted log file
state_content = mytextread(fullfile(state_path, statename));

%read info from logfile
temp_log_info = Words2cell(log_content{2}); % separate the 2 words making up the 2nd line of log file
ExpID = temp_log_info{1}; 
identifier = temp_log_info{2};
%---------------------------------
%new logfiletypes(from april on)
%---------------------------------
%old logfiletypes(before april); then you have to set ExpID here by hand
%ExpID = 'RG09097';
%identifier = char(temp_log_info{1}(1));

%read info from statefile
temp1= textscan(state_content{2},'%s','delimiter',':');%quality
temp2= textscan(state_content{6},'%s','delimiter',':');%adjusted unit
temp3= textscan(state_content{4},'%s','delimiter',':');%gain
rec_quality = char(temp1{1}(2));
rec_base_unit = char(temp2{1}(2));
rec_gain = str2num(char(temp3{1}(2)));

%open SGSR dataset and read info from there
try,
    infoSGSRset = sgsrdataset(ExpID,identifier);
    stim_property_unit = infoSGSRset.xunit;
    stim_property_values = infoSGSRset.xval;
catch,
    warning(lasterr);
    stim_property_unit = '???';
    stim_property_values = [];
end

%generate output struct
date = datestr(min([abffiles.datenum]),20); % date of oldest abf in the neighborhood
divider = '__________';
start_rec_str = datestr(start_rec,'HH:MM:SS');
stop_rec_str = datestr(stop_rec,'HH:MM:SS');
recording_infos = struct('Date_of_recording', date, 'Timestamp_first_abf', ...
    start_rec_str, 'Timestamp_last_abf',stop_rec_str,'d_____',divider,'SGSRExpID',...
    ExpID,'SGSRidentifier',identifier, 'Stim_Property_Unit',stim_property_unit, ...
    'Stim_Property_Values',stim_property_values,'Rec_Quality',rec_quality, ...
    'Rec_Base_unit',rec_base_unit,'Rec_Gain',rec_gain);

















