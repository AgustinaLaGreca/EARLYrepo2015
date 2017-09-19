function JSlinkaxes(ah, axdim);
% JSlinkaxes - link axes using common sense
if ismember('x', lower(axdim)),
    XL = [inf -inf];
    for ii=1:numel(ah),
        axes(ah(ii));
        XL = [min(XL(1), min(xlim)), max(XL(2), max(xlim))];
    end
end
if ismember('y', lower(axdim)),
    YL = [inf -inf];
    for ii=1:numel(ah),
        axes(ah(ii));
        YL = [min(YL(1), min(ylim)), max(YL(2), max(ylim))];
    end
end


linkaxes(ah, axdim);

if ismember('x', lower(axdim)),
    xlim(XL);
end
if ismember('y', lower(axdim)),
    xlim(YL);
end

figh = parentfigh(ah(1));
POS = get(figh, 'position');
set(figh, 'position', POS*0.999); % force vivibility of legends ..



