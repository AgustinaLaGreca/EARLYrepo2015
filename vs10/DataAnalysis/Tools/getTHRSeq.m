function [seqNr, dsTHR] = getTHRSeq(Exp, cellNr)
% getTHRSeq - Get THR sequence for this particular data file and cell number
%
% [seqNr, dsTHR] = getTHRSeq(dataFile, cellNr)
%
% Get THR sequence for this particular data file and cell number, without
% using the user data. The THR dataset object is returned as second output parameter.

%B. Van de Sande 29-08-2003
%K. Spiritus

% ---------------- CHANGELOG -----------------------
%  Mon Jan 17 2011  Abel	* Check if the found sequence is really a THR
%								curve by looking at the stimtype.
%							* The THR dataset is now returned als output.
%							* Cleaned code / rewrite function
%  Fri Jan 21 2011  Abel   
%   - bugfix, when no THR found, set dsTHR to NaN

%% ---------------- Default parameters --------------
force = 1;			% Force caching to off

%% ---------------- Main program --------------------
% Get all THR recordings from this cell
thr_recordings = local_get_thr(Exp,cellNr);

% Return [] if nothing found 
if isempty(thr_recordings)
	seqNr = NaN;
	dsTHR = NaN;
 	warning('EARLY:Debug', 'Unable to get THR sequence for cell:%d', cellNr);
	return;
end

% Get the latest THR recorded == CONVENSION
index_THR = -1;
highestRecFromCell = -1;
for i=1:length(thr_recordings)
   RecFromCell = elem(strsplit(thr_recordings{i}.dsID,'-'),2);
   if str2num(RecFromCell)  > highestRecFromCell
       highestRecFromCell = str2num(RecFromCell);
       index_THR = i;
   end
end

if highestRecFromCell == -1
   error('EARLY:ERROR','Unable to find THR for cell number:%d',cellNr); 
end

% Load the dataset for the THR data
dsTHR = read(dataset,Exp,str2num(thr_recordings{i}.ds_unique_number));
seqNr = str2num(thr_recordings{i}.ds_unique_number);


end
%% ---------------- Local functions -----------------
%{
	get_max_thr_sequence_nr_:
	Find last sequence number in list which contains 'th' in the
	indetifier
%}
function idx = get_thr_sequence_nr_(idList, allCellNrs, cellNr)
	idx = intersect(find(allCellNrs == cellNr), strfindcell(lower(idList), 'th'));
end

function thr_recordings = local_get_thr(Exp,cellNr)
    thr_recordings = {};
    % Open the log file
    log_file_name = ['C:\ExpData\Exp' filesep Exp filesep Exp '.log'];
%    log_file_name = [ 'C:\ExpData\Exp' filesep Exp.ID.Name filesep Exp.ID.Name '.log'];
    fid = fopen (log_file_name,'rt');

    if fid == -1
       error('EARLY:Error','Unable to open the log file for experiment: %s',Exp); 
    end
    
    % Read line by line and check if it has info on THR
    thr_lines = {};
    line = fgetl(fid);
    while ischar(line)
        if strfind(line,'THR')
           if strfind(line,'#') % Dataset declarations contain a #
              thr_lines{end+1} = line; 
           end
        end
        line = fgetl(fid);
    end
    
    % Parse the line and check if it is from this cell
    for i=1:length(thr_lines)
       dsID = elem(strsplit(thr_lines{i},'  '),2);
       CellNr = elem(strsplit(dsID,'-'),1);
       if str2num(CellNr) == cellNr
          ds_unique_number = elem(strsplit(elem(strsplit(elem(strsplit(thr_lines{i},'  '),1),' '),3),':'),1); 
          thr_struct = CollectInStruct(dsID,ds_unique_number);
          thr_recordings{end+1} = thr_struct;
       end
    end
    
    
end
