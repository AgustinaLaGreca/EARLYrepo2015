function D=dummy(VD, ExpName, ids, Stim);
% dataset/dummy - dummy dataset having given experiment name and rec index
%    D=dummy(dataset(), 'RG11303', 23) is a dummy dataset for which
%        expname(D) equals 'RG11303'
%        and irec(D) equals 23.
%
%    dummy(dataset(), ExpName, I) returns an array of summy datasets if I
%    is an array.
%
%    dummy(DS) is the same as dummy(dataset(), expname(DS), irec(DS))
%
%    D=dummy(dataset(), ExpName, I, Stim) also provides D with a stimulus
%    field Stim.
%
%    See also dataset/isdummy, experiment/dummy.

Stim = arginDefaults('Stim', struct());

if nargin==1,
    for ii=1:numel(VD),
        D(ii) = dummy(dataset(), expname(VD(ii)), irec(VD(ii)));
    end
    D = reshape(D, size(VD));
    return;
end

if numel(ids)>1,
    for ii=1:numel(ids),
        D(ii) = dummy(VD, ExpName, ids(ii));
    end
    D = reshape(D, size(ids));
    return;
end

% ----single irec from here----
D = dataset(); % void
D.ID.Experiment = dummy(experiment(), ExpName);
D.ID.iDataset = ids;
D.ID.uniqueID = 0; % indicating this is a dummy
D.ID.iCell = nan;
D.ID.iRecOfCell = nan;
D.Stim = Stim;


