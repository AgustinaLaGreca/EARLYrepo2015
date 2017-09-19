function X=xspikeAlign(X, T0, T1, dtMax);
% xspikeAlign - time-align individual events with their template
%    X=xspikeAlign(X, t0, t1, dtMax), where X is the output of xspikes, aligns 
%    the indivual events of category K to best match the segment of their
%    template between t0(K) and t1(K) ms. 
%
%    See also xspikes, xspikeplot.


if nargin<2, T0=-inf; end
if nargin<3, T1=inf; end
if nargin<4, dtMax = 0.4; end; % ms max time shift

PrecDur = 0.5; % ms default interval preceeding max(template) to be used for alignment

[T0, T1] = SameSize(T0, T1, ones(1,X.Ncat));
dt = X.dt;
[Tev0, Tev1] = minmax(X.t0snip);
% get sample range of template segment
I0 = 1+round((T0-Tev0)/dt); % first sample of segment
NsamSeg = round((T1-T0)/dt);
dImax = round(dtMax/dt);
I1 = I0 + NsamSeg-1; % last sample of segment
for icat=1:X.Ncat,
    ievent = find(X.catIdx==icat);
    Ev = getevent(X,ievent);
    NsamEv = size(Ev,1);
    if isempty(Ev), continue; end
    Template = median(Ev,2);
    if isinf(T0(icat)), % segment is PrecDur ms preceeding max of template
        [dum, i1] = max(Template);
        i0 = i1+1-round(PrecDur/dt);
    else, % cut out segment as specifiied by caller
        i0 = I0(icat); i1 = I1(icat);
    end
    Seg = Template(i0:i1); % specified segment of tempolate to be used for for alignment
    i0ev = max(1,i0-dImax(1));
    i1ev = min(NsamEv, i1+dImax(end));
    Ev = Ev(i0ev:i1ev,:);
    Ishift = locateSegment(Seg,Ev)+i0ev-i0;
    tshift = dt*Ishift;
    X.tmax(ievent) = X.tmax(ievent)+tshift';
end














