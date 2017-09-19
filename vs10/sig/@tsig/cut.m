function T = cut(T, t0,t1,tstart,defVal);
% tsig/cut - cut interval from tsig object
%
%   R=T.cut(t0,t1) or R=cut(T,t0,t1) cuts an interval from T starting at t0
%   ms and ending at t1 ms. More accurately, the interval starts at t0 ms
%   and its duration is (t1-t0) ms. By convention, R.starttime is 0 ms.
%   Note that the temporal position of the cut within the waveform of T 
%   depends on T.onset. If t0 and/or t1 our outside the interval 
%   (T.starttime .. T.endtime) then zeros are padded to obtain at the 
%   requested duration of R. This convention guarantees the duration of
%   the cut to be independent of the time domain of T. Note that the cut 
%   can only contain an integer number of cycles, so t0 and t1 are rounded
%   toward the nearest multiple of T.dt.
%
%   Special values for t0 and t1 are:
%     t0 = -inf  means T.starttime 
%     t1 =  inf  means T.endtime
%
%   R=cut(T,t0,t1,tstart) or R=T.cut(t0,t1,tstart) sets the onset of R to
%   tstart [ms]. The default values is 0.
%
%   t0 and/or t1 and/or t3 may be vectors, in which case their elements are 
%   applied to the respective channels of T. Note that specifying multiple-valued
%   t0,t1, or tstart to a single-channel T results in a multi-channel R.
%
%   R=cut(T,t0,t1,tstart,defVal) or R=T.cut(t0,t1,tstart,defVal) pads 
%   defVal instead of zeros. defVal must be a numerical or logical scalar.
%
%   T.cut(t0,t1)=B  or  T=cut(T,t0,t1,cutting(B)) is a cut-and replace call, 
%   which returns T with the segment T.cut(t0,t1) replaced by B. For this
%   cut-and-replace call t0 and t1 may not exceed the (channelwise) time range of T.
%   If B is a numeric or logical value, the segment T.cut(t0,t1) is filled
%   with that constant value. B may also be a row vector whose length is
%   compatible with T.nchan, in which case each T(ichan).cut(t0,t1) is
%   are replaced by the constant B(ichan). If T is single-channel and B a
%   vector, a multi-channel tsig is created. Finally, B may also be a tsig
%   object, which must have the proper (per-channel) duration(s) t1-t0, and
%   a compatible number of channels. Again, a single-channel T is turned 
%   into a multi-channel tsig when B.nchan>1.
%   
%   See tsig/dur, tsig/starttime, tsig/endtime, tsig/dt, tsig/randcut.

if nargin<2, % probably from a subsref construction T.cut(..). Delegate to cutting/subsref to find the (..)
    T = cutting(T); % recursion in subsref should fetch the params t0, t1.
    return;
end 
if nargin<3, 
    if isTsig(t0),
        % probably from a subsasgn construction T.cut(..)=X. The real work has ...
        % ... already been done by recursion in tsig/subsasgn recursion. Return ...
        % ... the second arg because it contains T with the segment replaced
        T = t0;
        return;
    end
end
if nargin<4,
    tstart=0; % default onset of the output is 0 ms.
end
if nargin<5,
    defVal=0; % pad zeros in out-of-range time intervals.
else,
    if ~isscalar(defVal) || ~isNumOrLogical(defVal),
        error('defVal argument must be logical or numerical scalar.');
    end
end

% check t0 and t1 args, and expand their special values.
error(dimensionTest(t0,'row', 'Input argument t0'));
error(dimensionTest(t1,'row', 'Input argument t1'));
error(numericTest(t0,'real/nonnan', 'Input argument t0 is '));
error(numericTest(t1,'real/nonnan', 'Input argument t1 is '));
% t0 and t1 must have as many values as there are channels in T
Nchan = nchan(T);
try, [t0, t1] = SameSize(t0,t1,1:Nchan);
catch, error('Incompatible channel counts.')
end
% Expand special values. See help text.
Tonset = onset(T); Toffset = offset(T);
ihit = find(t0==-inf); t0(ihit) = Tonset(ihit); % -inf means T.onset
ihit = find(t1==inf); t1(ihit) = Toffset(ihit); % inf means T.offset
% At this stage, t0 and t1 should be finite numbers; check it.
error(numericTest(t0,'noninf', 'Input argument t0 is '));
error(numericTest(t1,'noninf', 'Input argument t1 is '));

% if the last arg is a cutting object, by definition we're dealing with a
% cut-and-replace call. See help text.
if isa(tstart, 'cutting'),
    Rep = cuttee(tstart); % unwrap the replacer
    [T, Mess] = local_replace(T, t0,t1, Rep);
    error(Mess);
    return;
end
% if we're still here, we're doing a "passive" cut: no replacement.
% tstart is the onset of the result. Test it.
try, tstart = SameSize(tstart, 1:Nchan);
catch, error('Incompatible channel counts.')
end
error(dimensionTest(tstart,'row', 'Input argument tstart'));
error(numericTest(tstart,'real/nonnan/noninf', 'Input argument tstart is '));

% get the sample indices
dT = dt(T); % sample period in ms
t2i = @(t,ichan)round(t(ichan)/dT); % portable time-to-index converter
Nsam = nsam(T); % #samples in each channel
for ichan=1:Nchan,
    n0 = 1+t2i(t0-Tonset,ichan); % starting index in waveform. 1+  because Matlab base-1 convention
    N = t2i(t1-t0,ichan); % sample count
    n1 = n0+N-1; % ending sample in waveform
    W = repmat(defVal,N,1); % out-of-range values are zero
    i0 = max(1,n0);
    i1 = min(Nsam(ichan), n1);
    Nsubst = i1-i0+1; % #samples to be substituted
    W(i0-n0+(1:Nsubst)) = T.Waveform{ichan}(i0:i1);
    T.Waveform{ichan} = W;
    T.t0(ichan) = tstart(ichan);
end

%---------------------------------------------------

function [T, Mess] = local_replace(T, t0,t1,Rep);
% cut-and-replace. Replace segment T.cut(t0,t1) by Rep.
Mess = ''; % optimistic default
if isempty(Rep), % the segment will be removed, not replaced
    Nchan_rep = 1; % arbitrary value to suppress incompatible #chan error further down
elseif isNumOrLogical(Rep), % rep must be scalar or row vector
    if ~isequal(1,size(Rep,1)),
        Mess = 'In an assignment T.cut(t0,t1)=X, any numeric or logical X must be a row vector, scalar, or [].';
        return;
    end
    Nchan_rep = size(Rep,2);
elseif isTsig(Rep),
    Nchan_rep = nchan(Rep);
else,
    Mess = 'RHS of assigment  T.cut(t0,t1)=X  must be numeric, logical, or a tsig object';
    return;
end
Nchan = [nchan(T) Nchan_rep];
if length(unique([1 Nchan]))>2,
    Mess = 'Incompatible channel count in cut-and-replace assigment T.cut(t0,t1)=X.';
    return;
end
Nchan = max(Nchan);
% if this call increases the #channels of T (see help text), duplicate the channels of T 
if nchan(T)<Nchan, T = repmat(T,1,Nchan); end
% now for the timing
[t0,t1] = SameSize(t0,t1,1:Nchan);
% convert t0 and t1 to indices and test for exceeding T's range
dT = dt(T); % sample period in ms
t2i = @(t)round(t/dT); % portable time-to-index converter
i0=1+t2i(t0-onset(T));
i1=t2i(t1-onset(T));
if any(i0<0) || any(i1>nsam(T)),
    Mess = 'Selected time interval in cut-and-replace call T.cut(t0,t1)=X exceeds limits of T.'
    return;
end
nsamrep = i1-i0+1; % channelwise sample count of segment-to-be-replaced
nsamrep = max(0,nsamrep); % # samples cannot be negative
% now we're ready for the actual substitution
ir = [i0;i1]; % LHS index range in 2xNchan matrix
if isempty(Rep), % remove the segment; no checking needed here
    for ichan=1:Nchan, % remove the segment
        T.Waveform{ichan}(vcolon(ir(:,ichan))) = [];
    end
elseif isNumOrLogical(Rep),
    for ichan=1:Nchan,
        T.Waveform{ichan}(vcolon(ir(:,ichan))) = Rep(min(ichan,end));
    end
else, % Rep is a tsig
    if any(nsamrep~=nsam(Rep)),
        Mess = 'In an assigment T.cut(t0,t1)=S, any tsig object S must have the same # samples as T.cut(t0,t1).'
        return
    end
    for ichan=1:Nchan,
        if nsamrep(ichan)>0
            T.Waveform{ichan}(vcolon(ir(:,ichan))) ...
                = Rep.Waveform{min(ichan,end)};
        end;
    end
end




