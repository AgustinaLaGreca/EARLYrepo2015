function ds = fastcache(G, Irec, ds);
% getds/fastcache - fast caching of datasets accessed by getDS object
%   fastcache(G, irec, ds) stores ds in entry #irec of G's cache.
%   Arrays irec are allowed. 
%   
%   ds = fastcache(G, irec) retrieves dataset from #irec of G's cache.
%   Arrays irec are allowed. Entries not present in the cache wil result in
%   void dataset entries.
%
%   
%   See also getDS, getDS/subsref.

persistent CCC % the cache

if isempty(CCC) || ~isfield(CCC, G.XP),
    CCC.(G.XP) = emptystruct('irec', 'ds');
end

if nargin==3, % store one by one
    for ii=1:numel(Irec),
        irec = Irec(ii);
        if ~any([CCC.(G.XP).irec]==irec), % not stored yet
            CCC.(G.XP) = [CCC.(G.XP), CollectInStruct(irec, ds)];
            if numel(CCC.(G.XP))>G.Nmax, % remove the first (=oldest) entry
                CCC.(G.XP)(1) = [];
            end
        end
    end
else, % retrieve one by one
    for ii=1:numel(Irec),
        irec = Irec(ii);
        if any([CCC.(G.XP).irec]==irec), % return cached ds
            ds(ii) = CCC.(G.XP)([CCC.(G.XP).irec]==irec).ds;
        else, % return void dataset
            ds(ii) = dataset();
        end
    end
end






