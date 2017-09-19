function Y = repmat(varargin);
% playsig/repmat - produce error to prevent user from constructing playlist arrays

error('Invalid repmat. Playlist objects may not be placed in arrays or matrices. Use cell arrays instead.');


