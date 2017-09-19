function RueXlswrite(S,FN);
% RueXlsWrite - write Rue struct to XLS file
%    RueXlswrite(S,FN);

Nrow = numel(S);
Fnames = fieldnames(S);
Ncol = numel(Fnames);
xlswrite(FN , Fnames');
Alphabet = {}; 
for ii=1:26, Alphabet{ii} = char(ii+64); end
for ii=1:26, Alphabet{ii+26} = ['A' char(ii+64)]; end
for ii=1:26, Alphabet{ii+2*26} = ['B' char(ii+64)]; end
for icol=1:Ncol,
    fname = Fnames{icol};
    if islogical(S(1).(fname)), 
        [S.(fname)] = DealElements(double([S.(fname)]));
    end
    Values = {S.(fname)}';
    RangeStr = [Alphabet{icol} '2:' Alphabet{icol} num2str(Nrow+1)];
    xlswrite(FN, Values, RangeStr);
end



