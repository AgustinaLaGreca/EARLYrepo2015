function D = subsref(G, S);
% getds/subsref - make Getds object get a dataset
%   G(4) returns dataset # 4 of G's experiment.
%   G(3:5) and G([2 5 8]) return datset arrays
%
%   P = G({1 [2 6]}) returns a 1x2 pooled_dataset array equal to
%       [pooled_dataset(G(1)) pooled_dataset(G([2 6]))]

if length(S)>1, % use recursion from the left
    y = subsref(G,S(1));
    D = subsref(y,S(2:end));
    return;
end

% single-elem S from here
switch S.type,
    case '()',
        if numel(S.subs)>1,
            error('Multiple indices not allowed for Getds objects.');
        end
        idataset = S.subs{1};
        doPool = iscell(idataset);
        if doPool,
            D = [];
            for jj=1:numel(idataset),
                idcell = idataset{jj};
                d = fastcache(G, idcell);
                for ii=1:numel(d),
                    if isvoid(d(ii)),
                        d(ii) = read(dataset, G.XP, idcell(ii));
                        fastcache(G, idcell(ii), d(ii));
                    end
                end
                D = [D pooled_dataset(d)];
            end
        else,
            D = fastcache(G, idataset);
            for ii=1:numel(D),
                if isvoid(D(ii)),
                    D(ii) = read(dataset, G.XP, idataset(ii));
                    fastcache(G, idataset(ii), D(ii));
                end
            end
        end
    otherwise,
        error('Invalid use of Getds object. The only valid syntax is G(k).');
end





