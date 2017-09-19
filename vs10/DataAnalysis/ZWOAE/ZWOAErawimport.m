function W = ZWOAErawimport(expID, recID)
%ZWOAErawimport -import ZWOAE recording that hold raw (unprocessed) signal
%
%    W = getZWOAEdata(expID, recID) returns datastruct stored in
%    the raw datafile with extension 'ZWOAE' holding the data of recording
%    recID of gerbil expID.
%
%    For file names, see ZWOAEfilename.
%
% Multiple recID are allowed, the returned variable W is then a struct array
%
%
% See also ZWOAEimport, getZWOAEdata, ZWOAEfilename

%---------- handle multiple recID recursively ------------
if numel(recID) > 1,
    for ii = 1:numel(recID),
        W(ii) = ZWOAErawimport(expID, recID(ii));
    end
    return
end


%---------- single recording from here------------
[FN, DD] = ZWOAEfilename(expID, recID,[],'raw');
FFN = fullfile(DD,FN);
if ~exist(FFN,'file'),
    error(['Raw datafile ''' FN ''' not found in ZWOAE data path.']);
end
tmp = load(FFN,'-mat');
if isstruct(tmp), %data was saved as a struct already...
    fns = fieldnames(tmp); 
    W = tmp.(fns{1}); %strip one "struct-layer" of from import
else
    W = tmp; %data was not saved as a struct, no need to strip one layer off
end