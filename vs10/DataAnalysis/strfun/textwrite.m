function textwrite(fid, TXT);
% textwrite - write text matrix to open textfile; deblank each line of text

Nline = size(TXT,1);
for iline=1:Nline,
   lin = deblank(TXT(iline,:)); 
   %fprintf(fid, '%s\n', lin);  % commented by Hsin-Wei on 03.Oct.17
   %% following added by Hsin-Wei on 03.Oct.17
   fileID = fopen(fid,'a'); % convert fid from string to a file ID so that fprintf can write text to the file
   try
       fprintf(fileID,'%s\n', lin); % somehow lin sometimes becomes a CELL and causes error in fprintf
   catch err
       lin = cell2words(lin);       % convert lin from CELL to CHAR so that fprintf can work
       fprintf(fileID, '%s\n', lin);
   end
   fclose(fileID);
end
