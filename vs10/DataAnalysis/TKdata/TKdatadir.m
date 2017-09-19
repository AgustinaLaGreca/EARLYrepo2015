function DD = TKdatadir(Exp);
%  TKdatadir - datadir for TK's ABF data.
%
%  TKdatadir returns full path for datadir for TK's ABF data. The precise
%  path depends on the computer.
%
%  TKdatadir('JL') returns full path for datadir for JL's ABF data. 
%  TKdatadir('TK') is the same as TKdatadir().
%
% See also TKpool.

if nargin<1,
    Exp = 'TK';
end
[Exp, Mess] = keywordMatch(upper(Exp), {'TK' 'JL'}, 'hybrid experimenter');
error(Mess);

switch Exp
    case 'TK',
        switch CompuName,
            case 'LHAW',
                %DD = 'D:\USR\Marcel\sortTKdata\tempABFstore';
                DD = 'D:\USR\Thomas\Data';
            case {'SIUT', 'CLUST'},
                DD = 'D:\RawData\TKdata';
            case 'Cel',
                DD = 'C:\D_Drive\rawData\ThomasABFdata';
            case 'KULAK',
                DD = 'D:\USR\Thomas\data';
            case 'Ee1285a',
                DD = 'D:\rawdata\TKdata';
            otherwise,
                error('Cannot find TK data storage on this computer.');
        end
    case 'JL',
        switch CompuName,
            case 'LHAW',
                DD = 'D:\USR\Jeannette\Data';
            case {'SIUT', 'CLUST'}
                DD = 'D:\RawData\JLdata';
            case 'Cel',
                DD = 'C:\D_Drive\rawData\JeannetteABFdata';
            case 'KULAK',
                DD = 'D:\USR\Jeannette\data';
            case 'Ee1285a',
                DD = 'D:\rawdata\JLdata';
            otherwise,
                error('Cannot find JL data storage on this computer.');
        end
end



