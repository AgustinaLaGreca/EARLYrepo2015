function n=nchan(C);
% stimcontext/nchan - number of active DA channels from stimulus context
%     nchan(S) - returns the number of active DA channels from stimulus 
%     context S. 
%
%     See stimcontext.

switch lower(C.DAchan),
    case {'left' 'right'}, n=1;
    case {'both'}, n=2;
end
