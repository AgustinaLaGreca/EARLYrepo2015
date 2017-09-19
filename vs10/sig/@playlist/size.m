function S = size(P, idim)
% playlist/size - SIZE for playlist objects.
%
%   Tricky! Size(P) is by definition 1xNsamplay(P).
%
%   See playlist.

S = [nsamplay(P) 1];
if nargin>1,
    if idim<3, S=S(idim);
    else, S=1;
    end
end







