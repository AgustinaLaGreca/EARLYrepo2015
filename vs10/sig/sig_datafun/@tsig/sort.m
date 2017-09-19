function Y = sort(T, dim, Mode);
% tsig/sort - SORT for tsig objects.
%    Sort(S) for tsig object S sorts the samples in each channel of S in
%    ascending order. Note that the same (idiotic) rules are used for
%    sorting complex numbers as described in Matlab SORT.
%  
%    Sort(S,1) is the same as sort(S). Thus by convention, the "first 
%    dimension" of tsig objects is the time dimension. This is consistent 
%    with tsig/vertcat appending tsig objects in time.
%
%    T=sort(S,2) operates cross channels. A tsig object T is returned 
%    in which the samples of S at each time instant are redistributed over 
%    the channels such that at any time T(1)<T(2)<.. Importantly, the 
%    channels are time-aligned before being processed. 
%
%    Sort(S,dim,'descend') sorts in descending order. Sort(S,dim,'ascend')
%    is the same as Sort(S,dim).
% 
%    See also tsig/mean, tsig/meadian, tsig/plus, tsig/exist, tsig/isvoid.

if nargin<2, dim=1; end
if nargin<3, Mode='ascend'; end

if ~isequal(1,dim) && ~isequal(2,dim),
    error('Dim argument must be either 1 or 2.');
end
if isvoid(T),
    error('Invalid input: void tsig object.'); 
end

if dim==1, % sort each channel
    my_sort = @(x)sort(x,1,Mode); % 2nd arg = DIM
    Y = wavefun(my_sort, T, 'tsig');
elseif dim==2, % sort channels samplewise
    Nchan = nchan(T);
    E = exist(T); %#ok<EXIST> MLint does not realize exist is overloaded
    for epi=E, % paste channels living in this episode into matrix
        epi
        w = cut(T,epi.tstart, epi.tend);
        wm = [w.Waveform{epi.ichan}]; % waveform matrix
        wm = sort(wm,2,Mode); % sort along rows = channels
        for ii=1:length(epi.ichan), % insert sorted row into old position
            ichanLHS = epi.ichan(ii); % channel count in T
            ichanRHS = ii; % column count in wm matrix
            i0LHS = epi.i0(ii); i1LHS = epi.i1(ii);
            T.Waveform{ichanLHS}(i0LHS:i1LHS) = wm(:,ichanRHS);
        end
    end
    [Y, Mess] = test(T);
    error(Mess);
end



