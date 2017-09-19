function Y = DPampCoeff(N);
% DPampCoeff - analyze DP amplitude differences from combinatorics
%    DPampCoeff(N) displays a list of types of Nth-order terms and their
%    multinomial  coeffiecnts.
%
%    See also uberzwuis.

p = partitions(N); % rows are the exponents of the subsequent factors
% discard the identity of the factors, only keep the combinatoric
% differences
MM=unique(sort(partitions(N),2), 'rows');
Nconfig = size(MM,1); % # configurations, e.g., a^2bc^3
for irow=1:size(MM,1), 
    Row = MM(irow,:); 
    Coeff(irow) = factorial(sum(Row))/prod(factorial(Row));
    Config{irow} = Row(Row>0);
end
% Keep only those configs for which sumweight=1 can be realized
iok = [];
for ic=1:Nconfig,
    c = Config{ic};
    Nc=numel(c);
    q = c(1);
    for ii=2:Nc,
        [q dq] = SameSize(q,c(ii)*[-1 1]);
        q = q+dq;
        q = q(:);
    end
    iok(ic) = ismember(1, abs(q));
end
iok = logical(iok);
[Coeff, Config] = deal(Coeff(iok), Config(iok));
Coeff = Coeff/mgcd(Coeff);
[Config, Coeff] = sortAccord(Config, Coeff, Coeff);
Coeff_dB = A2dB(Coeff/min(Coeff));

% display if requested
Nconfig = numel(Config);
if nargout<1,
    ABC = char('a'+(0:25));
    str = '';
    for ii=1:Nconfig,
        c = Config{ii};
        str = [str num2str(Coeff(ii)) '*'];
        for jj=1:numel(c);
            str = [str ABC(jj)];
            if c(jj)>1,
                str = [str '^' num2str(c(jj))];
            end
        end
        str = [str ' + '];
    end
    disp(['  --------' nthstr(N) '-order DPs ---------'])
    str = str(1:end-3);
    disp(str)
    dBstr0 = [trimspace(num2str(Coeff_dB,2)) ' dB'];
    dBstr0 = strrep(dBstr0, ' ', '   ');
    dBstr = trimspace(num2str(diff(Coeff_dB),2));
    dBstr = ['+' strrep(dBstr, ' ', '   +')];
    disp(dBstr0);
    disp(['   '  dBstr ' dB']);
else,
    dB_Inc = diff(Coeff_dB);
    Y = CollectInStruct(Nconfig, Coeff, Coeff_dB, dB_Inc, Config);
end










