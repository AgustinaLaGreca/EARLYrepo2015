function Info = sys3CircuitInfo(Dev, InfoType);
% sys3CircuitInfo - info on loaded circuit 
%
%   sys3CircuitInfo(Dev) returns a struct with fields:
%
%            Device: Name of the device
%       CircuitFile: Filename of the loaded circuit
%              Fsam: sample rate o circuit in kHz
%        CycleUsage: cycle usage in %
%         Component: list of generic circuit components
%          ParTable: list of ParTable components
%           SrcFile: list of SrcFile components
%            ParTag: info on ParTags in struct array with fields
%                    Name, DataType, TagSize, TagValue.
%
%   Items = sys3CircuitItems(Dev, InfoType) returns only the requested type
%   of circuit info. InfoType is a character array (string) indicating any
%   of the fields listed above. These are case-insensitive and may be 
%   uniquely abbreviated (e.g. the abbreviation ParT is not unique).
%   Empty values are returned if no circuit is loaded to Dev.
%
%   Examples:
%
%     ItemStruct = sys3CircuitInfo;
%     ItemStruct = sys3CircuitInfo('RP2_1');
%     ItemStruct = sys3CircuitInfo('RP2_1','Comp');
%
%   Note: the info is not directly retrieved from the TDT devices, but 
%   from a set of book keeping variables generated by sys3loadCircuit.
%
%   See also sys3status, sys3Fsam, sys3ParTag, sys3loadCircuit.

if nargin<3, flag = ''; end % default: no TDT polling; ask circuitInfo

%default to sys3defaultdev if no Dev is specified
if nargin < 1, Dev = ''; end

%Check whether a circuit is loaded on specified device.
%If not, exit the function and return empty variable Items
%Use CircuitInfo for this test to avoid unnecessary communication with TDT
%devices.
if isempty(private_CircuitInfo(Dev)), 
    Info = []; 
    return
end

% get all info
Info = private_CircuitInfo(Dev);

% select specific info if requested
if nargin>1 && ~isempty(InfoType), % expand Type from List
    InfoTypeList = fieldnames(Info); % available types of info
    [InfoType, Mess] = keywordMatch(InfoType, InfoTypeList);
    error(Mess);
    Info = Info.(InfoType);
end






