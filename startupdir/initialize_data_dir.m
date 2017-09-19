function initialize_data_dir()
%initialize_data_dir sets the data directory for the user

% Get the location of the folder where the Early data is stored
folder_name = uigetdir('C:\','Choose the folder with all Early Experiments:');

if isempty(folder_name) || folder_name == 0
    error('No data folder selected! Please restart Early');
end

datadir(folder_name);

end

