function [T, Y, isort]=localmax(t,y, minsegdur);
% localmax - determine local maxima of array
%   [T, Y, isort]=localmax(t, y, minsegdur) determines the local maxima Y of y and
%   their time of occurence T. A local maximum of y is a value y(k) with
%   y(k)>y(k-1) && y(k)>y(k+1).
%
%   Inputs
%        t: time axis or, if single number, its spacing from zero start
%        y: trace to be analyzed (same size as t)
%     minsegdur: local maxima at tmax are only accepted if in the interval
%          (tmax-minsegdur/2, tmax+minsegdur/2) no sample exceed the value
%          at tmax. If minsegdur is unspecified or set to [], no such check
%          is performed.
%   Outputs
%        T: time where local maxima of y occur
%        Y: local maxima themselves
%    isort: index array sorting Y in descending order.
%
%   See also MinMax, DIFF.

if nargin<3, minsegdur=[]; end

if ~isvector(t) || ~isvector(y),
    error('Both t and y input args must be vectors.');
end
t = t(:); y = y(:);
if isequal(1,numel(t)) && (numel(y)>1), % t is really spacing
    t = timeaxis(y, t);
end
if ~isequal(size(t),size(y)),
    error('Input args t and y must have same number of elements.');
end
% local max: 1st derivative changes sign
dy = diff(y); % first derivative
imax = 1+find( (dy(1:end-1)>0) & (dy(2:end)<=0) ); % +1 accounts for loss of one sample per diff operation
Y = y(imax); T = t(imax);
% only do sorting if requested - may be time-consuming
if nargout>2,
    [dum, isort] = sort(Y,'descend');
end

if ~isempty(minsegdur),
    dt = diff(t(1:2)); % time spacing
    isam_halfDT = round(0.5*minsegdur/dt);
    iok = true(size(imax));
    for ii=1:numel(imax),
        isam_pk = imax(ii); % sample index of peak
        i0 = max(1,isam_pk-isam_halfDT);
        i1 = min(numel(y),isam_pk+isam_halfDT);
        iok(ii) = all(y(isam_pk)>=y(i0:i1));
    end
    T = T(iok);
    Y = Y(iok);
    if nargout>2,
        [dum, isort] = sort(Y,'descend');
    end
end




