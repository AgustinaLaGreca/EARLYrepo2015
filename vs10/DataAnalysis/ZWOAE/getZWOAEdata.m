function D = getZWOAEdata(igerbil, idata, flag)
% getZWOAEdata - read ZWOAE data from .ZWcom file
%    D = getZWOAEdata(igerbil, irec) returns datastruct stored in
%    the datafile with extension 'ZWcom' holding the data of recording irec
%    of gerbil igerbil.
%
%    For file names, see ZWOAEfilename.
%
%    getZWOAEdata('RG07047_zDPs', [1 2 4]) returns datasets #1, 2 and 4 in
%    a struct array.
%
%    getZWOAEdata(igerbil, irec '-nosig') returns dataset(s) irec
%    with the signals stripped off, but all parameters in place. Such
%    dataless data can serve bookkeeping purposes.
%
%    See also ZWOAEdatadir, ZWOAEstruct.

if nargin<3, flag=''; 
elseif ~isempty(flag),
    [flag, Mess] = keywordMatch(flag, {'-nosig'}, 'flag argument');
    error(Mess);
end

if numel(idata)>1, % multiple datasets; handle recursively
    for ii = 1:numel(idata),
        D(ii) = getZWOAEdata(igerbil, idata(ii), flag);
    end
    return
end


%---------- single recording from here------------
[FN, DD] = ZWOAEfilename(igerbil, idata);
FFN = fullfile(DD,FN);
if ~exist(FFN,'file'),
    error(['Datafile ''' FN ''' not found in ZWOAE data path.']);
end
tmp = load(FFN,'-mat');
if isstruct(tmp), %data was saved as a struct already...
    fns = fieldnames(tmp); 
    D = tmp.(fns{1}); %strip one "struct-layer" of from import
else
    D = tmp; %data was not saved as a struct, no need to strip one layer off
end
if isequal('-nosig', flag), % strip off data
    D = rmfield(D, 'signal');
end
    
