function backup(Exp)
%backup will create a backup of the data on the server and add experiment and dataset info to the kquest database
%   Experiment is zipped and copied to the correct folder on the server. If
%   this was not succesfull the user will be informed. Next all datasets
%   are parsed and all experiment and dataset info is added to the database
%
%   Exp is the Exp which has to be backuped. Type= experiment object

%% Check if Experiment exists and contains datafiles
if ~isa(Exp,'experiment')
    error('Backup.m :: Input of the backup function is not an experiment object!');
end
exp_name = Exp.ID.Name;

if exist(experiment(), exp_name) == 0
    errordlg('There is no experiment currently loaded on the dashboard', ...
        'Error Non-existant Experiment');
    return;
end

exp_folder = folder(Exp);
exp_files = dir(exp_folder);


log_files = {};
expdef_files = {};
statuslog_files = {};
earlyds_files = {};
thr_files = {};
bin_files = {};

% remove the file and folder '.', '..', '...'
for i=1:length(exp_files)
   file_name = exp_files(i).name;
   if ~(strcmp(file_name,'.') || strcmp(file_name,'..') || strcmp(file_name,'...'))
       % Get the file extension
       file_extension = strsplit(file_name,'.');
       file_extension = file_extension{end};
       
       % Check if it is a log file
       if strcmpi(file_extension,'log')
           log_files{end+1} = file_name;
       end
       
       % Check if it is an ExpDef file
       if strcmpi(file_extension,'expdef')
           expdef_files{end+1} = file_name;
       end
       
       % Check if it is a Statuslog file
       if strcmpi(file_extension,'statuslog')
           statuslog_files{end+1} = file_name;
       end
       
       % Check if it is an EarlyDS file
       if strcmpi(file_extension,'earlyds')
           earlyds_files{end+1} = file_name;
       end
       
       % Check if it is a bin file
       if strcmpi(file_extension,'bin')
           bin_files{end+1} = file_name;
       end  
       
       % Check if it is a bin file
       if strcmpi(file_extension,'mat')
           thr_files{end+1} = file_name;
       end  
   end
end

if isempty(log_files) 
    errordlg('This experiment does not contain a log file. The backup could not be created.','Error no log file');
    return;
elseif isempty(expdef_files)
    errordlg('This experiment does not contain a ExpDef file. The backup could not be created.','Error no ExpDef file');
    return;
elseif isempty(statuslog_files)
    errordlg('This experiment does not contain a StatusLog file. The backup could not be created.','Error no StatusLog file');
    return;
end
%% Find the destination folder
user_name = username();

% DEBUG with local folder
SERVER_FOLDER = 'L:\early_data\Exp';
% SERVER_FOLDER = 'C:\usr\adriaan\server';

server_found = 0;
try
    cd(SERVER_FOLDER);
    server_found = 1;
catch
    errordlg('Unable to find the folder on the server');
    return;
end

%% Create a zip file
cd(exp_folder)
warning('off', 'MATLAB:zip:archiveName');
results = zip(exp_name,exp_folder);
warning('on', 'MATLAB:zip:archiveName');
if isempty(results)
   error('Early:EmptyZIP','Backup.m :: The zip file created contains no files'); 
   return;
end


%% Copy to server
if server_found
    [status,message,messageid] = movefile([pwd filesep exp_name '.zip'], SERVER_FOLDER);
end

% Check if file was moved
if status == 0
    errordlg('zip file was created but could not be moved to the server!','Can''t move zip file');
    error(messageid,message);
end
% Change back to the Early Sandbox folder
uiwait(helpdlg('Backup procedure was succesful!','Backup'));
cd('C:\Early\Sandbox')
%% Parse all data

% read each DS file one by one and collect them in a list of datasets
% ds_list = {};
% for i=1:length(earlyds_files)
%    
%     % get the dataset ID
%     ds_file = earlyds_files(i);
%     
%     splits = strsplit(ds_file{1},'_DS');
%     splits = strsplit(splits{2},'.');
%     
%     ds_id = str2num(splits{1});
%     
%     ds_list{end+1} = read(dataset(),exp_name, ds_id);
% end


%% Add data to the kquest database

%% Check if data was succesfully added to kquest

%% Let the user know backup was done succesfully


end

