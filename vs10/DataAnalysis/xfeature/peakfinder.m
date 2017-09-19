function [T, Y, isort, RM] = peakfinder(dt, y, minsegdur, MinVal);
% peakfinder - determine local maxima of array
%   [T, Y, isort] = peakfinder(dt, y, minsegdur) 
%   determines the local maxima Y of y and their times of occurrence T. 
%   A local maximum of y is a value y(k) with
%       y(k)>y(k-1) && y(k)>y(k+1).
%
%   Inputs
%        t: sample period of y 
%        y: trace to be analyzed (same size as t if t is an array)
%     minsegdur: local maxima at tmax are only accepted if in the interval
%          (tmax-minsegdur/2, tmax+minsegdur/2) no sample exceed the value
%          at tmax. If minsegdur<0, its absolute value is used, and among
%          any equal maxima occurring within this interval, the first is
%          selected.
%
%   Outputs
%        T: time where local maxima of y occur; same units as dt.
%        Y: local maxima themselves
%    isort: index array sorting Y in descending order.
%
%   [T, Y, isort] = peakfinder(dt, y, minsegdur, MinPeakVal) 
%   only returns those maxima of at least MinPeakVal.
%
%   See also Runmax.

[MinVal] = arginDefaults('MinVal', -inf);
DoCheck = minsegdur<0;
minsegdur = abs(minsegdur);
isrow = size(y,2)>1;
y = y(:);

NsamSeg = round(minsegdur/dt);
NsamSeg = max(3,NsamSeg);

[RM, imax] = runmax(y,NsamSeg);

qok = RM(imax)>=MinVal; % select largest peaks
imax = imax(qok);

if DoCheck,
    iclose = 1+find(diff(imax)<NsamSeg);
    while ~isempty(iclose),
        imax(iclose)=[];
        iclose = 1+find(diff(imax)<NsamSeg);
    end
end

T = dt*(imax-1);
Y = y(imax);
if nargout>2,
    [dum, isort] = sort(Y(:), 1, 'descend');
else,
    isort = nan;
end

if isrow,
    T = T.';
    Y = Y.';
    isort = isort.';
    RM = RM.';
end






