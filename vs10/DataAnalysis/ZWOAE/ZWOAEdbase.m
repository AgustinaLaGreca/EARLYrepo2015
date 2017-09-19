function D = ZWOAEdbase(igerbil,ilast,doStrip);
% ZWOAEdbase - display or generate databse of ZWOA parameters
%   ZWOAEdbase(31) displays a table containing the parameters of all
%   the recordings of gerbil # 31.
%
%   S = ZWOAEdbase(31) returns this info in a struct array and suppresses 
%   the screen display.
%
%  See also ZWOAEimport.

if nargin<2, ilast=nan; end
if nargin<3, doStrip=1; end

% find index of last recording, and use it unless ilast is prescribed
[FF DD ilastrecorded]=ZWOAEfilename(igerbil,inf);
if isnan(ilast), ilast = ilastrecorded; end

% get cached database
CFN = fullfile(DD, [mfilename dec2base(igerbil,10,3)]); Cacheparam = igerbil; Ncache=-1;
D = FromCacheFile(CFN, Cacheparam);
% find out how far we got last time when collecting the database
if isempty(D), imax = 0;
else, imax = max([D.recID]);
end
% anything new requested?
Nnew = ilast-imax; % # recordings to added 
%disp('-------')
%ilast, imax, Nnew
if Nnew>0,
    % supply newer data chunk by chunck to enable intermediate caching
    if Nnew<11,
        for ii=imax+1:ilast,
            D = [D, local_adrec(igerbil, ii)];
        end
        % store in cache
        ToCacheFile(CFN, Ncache, Cacheparam, D);
    else, % recursive call; add 10 at a time
        Nnext = imax;
        while Nnext<ilast,
            Nnext = Nnext+10
            D = ZWOAEdbase(igerbil,Nnext,0); % 0: don't strip
        end
    end
end

% strip off irrelevant fields and multiple-component field
if ~isempty(D) && doStrip,
    D = rmfield(D,{'Fzwuis' 'PHzwuis' 'PHsingle' 'PHsup' ...
        'a___________________' 'b___________________' 'micgain' 'dBV' 'dBSPL'});
end

% either display D or return it
if nargout<1, 
    structdisp(D, 'indexrow','off');
    clear D;
end

%====================================================
function d=local_adrec(igerbil, irec);
[FF DD]=ZWOAEfilename(igerbil,irec);
if ~exist(fullfile(DD,FF),'file'),
    d = [];
else,
    d = ZWOAEimport(igerbil,irec,'-nosig');
end




