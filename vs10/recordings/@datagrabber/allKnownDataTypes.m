function [L, gT] = allKnownDataTypes(G, Chan);
% datagrabber/allKnownDataTypes - list of all known data types
%   L = allKnownDataTypes(G) returns a cell string L whose elements are the 
%   names of known data types or "source devices", from whi dta are taken.
%
%   L = allKnownDataTypes(G, Chan) only returns those datatypes which may
%   be connected to a given channel Chan of the recording device. Chan is a
%   char string; currently available channels are
%      RX6_analog_1, RX6_analog_2, RX6_digital_1, PC_COM1.
%   
%   [L, grabType] = allKnownDataTypes(...) also returns a cell array
%   grabType whose k-th element is th class name of the datagrabber
%   subclass doing the grabbing for L{k}.
%   
%
%   See also datagrabber/dataType, datagrabber/getsourcedevicesettings.

Chan = arginDefaults('Chan', 'all'); % default: no selection
[Chan, Mess] = keywordMatch(Chan, {'all' 'RX6_analog_1' 'RX6_analog_2' 'RX6_digital_1' 'PC_COM1'});
error(Mess);

EE = {'grabevents'};
AA = {'grab_adc'};
CC = {'actabit'};

Ld1 = {'Spikes'};
La1 = {'Microphone-1' 'Microphone-2' 'Patch_pipette_1' 'Patch_pipette_2'  'RW_electrode' 'Micropipette' 'LaserDispl' 'LaserVelocity' 'DA-test'};
La2 = La1;
Com1 = {'Laser reflection'};

gTall = [repmat(EE, 1,numel(Ld1)) ...
    repmat(AA, 1,numel([La1 La2])) ...
    repmat(CC, 1,numel(Com1)) ];

switch Chan,
    case 'all',
        [L, gT] = deal([Ld1 La1 La2 Com1], gTall);
        [L iu] = unique(L);
        gT = gT(iu);
    case 'RX6_digital_1',
        [L, gT] = deal(Ld1, repmat(EE,1,numel(Ld1)));
    case 'RX6_analog_1',
        [L, gT] = deal(La1, repmat(AA,1,numel(La1)));
    case 'RX6_analog_2',
        [L, gT] = deal(La2, repmat(AA,1,numel(La2)));
    case 'PC_COM1',
        [L, gT] = deal(Com1, repmat(CC,1,numel(Com1)));
end









