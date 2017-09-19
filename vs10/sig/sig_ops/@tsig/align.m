function varargout=align(varargin);
% tsig/align - equalize dimensions and channel counts of tsig objects.
%    [S,T,Mess]=align(S,T) equalizes the time domains of S and T by padding zeros
%    outside their overlap region. Onsets of S and T are accounted for.
%    S and T must have compatible channels counts; if needed, channels will 
%    be duplicated so as to equalize the number channels in S and T.
%    Note that equalization is done in a channelwise way: S(1) is 
%    aligned with T(1), etc. Any differences across channels is not lifted.
%    If S and T are not suitable for alignment, a non-empty character
%    string Mess is returned, which is one of
%                 void: a void tsig was passed to binop
%                nchan: incompatible channel counts between S and T
%                 fsam: S and T have different sample frequencies.
%
%    [S,T,Mess]=align(S,T,defS,defT) fills the newly added portions of S with
%    the value defS instead of zeros; end the newly added portions of T
%    with the value defT. Both defS and defT must be scalar numeric or
%    logical values.
%
%    [S,T,..Z,Mess] = align(S,T,..Z, defS,...,defZ) aligns tsig objects S,T,..Z
%    while using respective default values defS,defT..defZ. Zeros are used
%    for missing def values.
%
%    See also tsig/arealigned.

Mess = '';
% find out what the args are about
Nsig = 0;
for ii=1:nargin,
    if isTsig(varargin{ii}), 
        Nsig=Nsig+1;
    else,
        break;
    end
end
% check suitability of tsig inputs for alignment. If something's wrong, return Mess & exit.
if any(cellfun(@isvoid, varargin(1:Nsig))),
    Mess = 'void';
end
Fsam = cellfun(@fsam,varargin(1:Nsig));
if length(unique(Fsam))>1,
    Mess = 'fsam';
end
% idem channel count
Nchan = cellfun(@nchan,varargin(1:Nsig));
if length(unique([1 Nchan]))>2,
    Mess = 'nchan'; 
end
Nchan = max(Nchan); % # channels of all output tsigs
varargout{Nsig+1} = Mess;
if ~isempty(Mess),
    return;
end
% provide zeros for missing def values 
if ~all(cellfun(@isNumOrLogical, varargin(Nsig+1:end))),
    error('Arguments must be tsig objects, optionally followed by their respective defValues.');
end
if ~all(cellfun(@isscalar, varargin(Nsig+1:end))),
    error('defValues must be scalar numeric or logical values.');
end
Ndef = nargin-Nsig;
defValues = num2cell(zeros(1,Nsig));
defValues(1:Ndef) = varargin(Nsig+(1:Ndef));
% At this stage, the signals are in varargin{1:Nsig} and the default values
% are in defValues{1:Nsig}

% deal with trivial cases
if Nsig==0,
    varargout={};
elseif Nsig==1, % nothing to align; simply return what's been passed
    varargout = varargin;
end

% ---the work starts here; most of it is bookkeeping---
% check fsam compatibility
% collect all onsets. Note that nchan may vary across inputs.
Onsets = cellfun(@onset,varargin(1:Nsig),'UniformOutput',false); 
[Onsets{:}] = SameSize(Onsets{:}); % duplicate single-channel ones if needed
Onsets = cat(1,Onsets{:}); % put in NsigxNchan matrix
% same story for offsets
Offsets = cellfun(@offset,varargin(1:Nsig),'UniformOutput',false); 
[Offsets{:}] = SameSize(Offsets{:}); % duplicate single-channel ones if needed
Offsets = cat(1,Offsets{:}); % put in NsigxNchan matrix
% the aligned signals run from min(onset) to max(offset) - channelswise (see help text)
t0 = min(Onsets);
t1 = max(Offsets);
% expand each tsig to have proper Nchan & time domain
for ii=1:Nsig,
    varargout{ii} = localExpand(varargin{ii},t0,t1,defValues{ii});
end

%---------------------
function S = localExpand(S,t0,t1,xdef);
Nchan = length(t0);
if nchan(S)<Nchan, % duplicate single channel
    S = repmat(S,1,Nchan);
end
t0S = onset(S);
t1S = offset(S);
DT = DT(S); 
nsamPre = round((t0S-t0)/DT); % # samples to be prepended, channelwise
nsamPost = round((t1-t1S)/DT); % # samples to be appended, channelwise
for ichan=1:Nchan,
    npre = nsamPre(ichan);
    npost = nsamPost(ichan);
    S.Waveform{ichan} = [repmat(xdef,npre,1); S.Waveform{ichan}; repmat(xdef,npost,1)];
end
S.t0 = t0;
% routine check if these changes are allowed
[S, Mess] = test(S);
error(Mess);








