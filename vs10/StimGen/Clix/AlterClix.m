function Y = AlterClix(flag, varargin);
%  AlterClix - alternating clicks a la Biedenbach & Freeman, 1964
%     AlterClix apens a GUI for presenting alternarting click stimuli

StartSPL = 30; % default dB start SPL
MinSPL = -20;
MaxSPL = 100;

flag = arginDefaults('flag', 'init');
if isequal('init', flag) && nargin<2,
    varargin = {StartSPL MinSPL MaxSPL};
end

Y = feval(['local_' flag], varargin{:});

%==============================
function figh = local_init(StartSPL, MinSPL, MaxSPL);
% GUI
TagName = 'AlterClixGUI';
figh = findobj('type', 'figure', 'tag', TagName);
if isempty(figh),
    CLR = [0.9130    0.6719    0.8968]; %1-0.4*rand(1,3);
    figh = figure;
    set(figh,'units', 'normalized', 'position', [0.419 0.692 0.238 0.118], 'tag', TagName, 'name', 'CLIX', ...
        'deletefcn', @local_stop, 'color', CLR, 'menubar', 'none', 'keypressfcn', @nope);
%     uicontrol('style', 'pushbutton', 'position', [ 50 50  80  30], 'string', 'PLAY', 'fontsize', 15, 'callback', @local_play);
%     uicontrol('style', 'pushbutton', 'position', [ 250 50  80  30], 'string', 'STOP', 'fontsize', 15, 'callback', @local_stop);
    uicontrol('style', 'text', 'position', [ 50 60  80  30], 'string', 'SPL', 'fontsize', 14, 'backgroundcolor', CLR);
    hspl = uicontrol('style', 'text', 'position', [ 220 60  120  30], 'string', [num2str(StartSPL) ' dB SPL'], 'fontsize', 14, ...
        'backgroundcolor', CLR, 'tag', 'SPL_reporter');
    uicontrol('style', 'slider', 'position', [130 60  80  30], 'min', MinSPL, 'max', MaxSPL, 'value', StartSPL, 'fontsize', 15, ...
        'callback', @(Src,Ev)local_spl(Src, hspl), 'backgroundcolor', CLR, 'tag', 'SPL_slider');
    uicontrol('style', 'text', 'position', [0 2  380  40], 'string', ...
        ['Connect (C1,C2) of PP16  to (CH1,CH2) of oscilloscope' char(10) 'and trigger on CH2 (=right ear = IPSI)'], 'fontsize', 11, 'backgroundcolor', 'w');
else,
    figure(figh);
end
% CIRCUIT
sys3loadCircuit('MSO_clicksearch', 'RX6');
sys3PA5(0);
h_sl = findobj(figh, 'tag', 'SPL_slider');
h_rp = findobj(figh, 'tag', 'SPL_reporter');
local_spl(h_sl, h_rp);
local_play;

function SPL = local_spl(Src, ht);
SPL =  get(Src, 'value');
Voltage = 9.9*dB2A(SPL-100);
set(ht, 'string', sprintf('%d dB SPL', round(SPL)));
sys3setpar(Voltage, 'ClickAmpl', 'RX6');

function local_play(Src, Ev);
sys3run('RX6');

function local_stop(Src, Ev);
sys3halt('RX6');





