function finish(flag)
% finish - logout function for EARLY
%   finish is called when quitting Matlab. It perfroms the following
%   actions:
%     - if TDT stuff is connected, set any attuators to 120 dB (this is
%       to prevent the crazy noise when turning of the computer, etc).
%    -  if TDT stuff is connected to RP or RX devices, clear any circuits.
%
%   See also startup, QUIT, sys3active, sys3loadCircuit.

try,
    if sys3isactive,
        disp('Resetting TDT devices ...')
        alldev = sys3devicelist;
        % PA5s - set to 120 dB
        iATT = strmatch('PA', alldev);
        for ii=1:length(iATT),
            ATT = alldev{iATT(ii)};
            invoke(sys3dev(ATT),'SetAtten',20);
        end
        % RP.. - clear any circuit
        iRP = strmatch('RP', alldev);
        for ii=1:length(iRP),
            RP = alldev{iRP(ii)};
            sys3loadCircuit('', RP);
        end
        % RX.. - clear any circuit
        iRX = strmatch('RX', alldev);
        for ii=1:length(iRX),
            RX = alldev{iRX(ii)};
            sys3loadCircuit('', RX);
        end
        disp('TDT devices reset.')
    end
catch,
    disp('Failing to reset TDT devices');
    pause(1);
end



