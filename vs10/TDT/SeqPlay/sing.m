function [SPinfo, Tone] = sing(Fsam, Dev);
% Sing - audio test using Seqplay functionality
%   Sing(Fsam) plays a 12-tone tune over each DAC channel using
%   a sample rate of ~Fsam kHz. Default sample rate is 10 kHz.
%   To repeat the song, type SeqplayGo.
%
%   sing(Fsam, Dev) uses TDT sys3 device Dev; default Dev is sys3defaultdev.
%
%   Sing('shuffle') realizes a new randomization and plays it.
%
%   See also Seqplay, SeqplayPlotStatus.

if nargin<2, Dev=''; end

if nargin<1, Fsam=10; end

% special case: reshuffle
if isequal('shuffle', lower(Fsam)),
    % re-specify play list
    Nmax = 6;
    if rand<0.5,
        seqplaylist([0 randperm(Nmax)], [3 2*ones(1,Nmax)], [randperm(Nmax) 0], [2*ones(1,Nmax) 3]);
    elseif rand<0.5,
        seqplaylist([],[], [randperm(Nmax) 0], [ones(1,Nmax) 3]);
    else,
        seqplaylist([randperm(Nmax) 0], [1+ones(1,Nmax+1)], [],[]);
    end
    %seqplaygo;
    return;
end

% initialize
%SPinfo = Seqplayinit(Fsam, Dev, '-divided');
%SPinfo = Seqplayinit(Fsam, Dev, '-triggering');
SPinfo = seqplayinit(Fsam, Dev);

% compute waveforms
RampDur = 5; % ms
Amp = 0.01;
Fsam = SPinfo.Fsam; % kHz
Dur = 200; % ms
Nsam = round(Fsam*Dur);
NsamRamp = round(Fsam*RampDur);
ramp = linspace(0,1,NsamRamp);
tt = linspace(0,200,Nsam);
for ii=1:12,
   freq = 0.5*2^(ii/12);
   Nt = round(Nsam/2*(1+rand));
   tone = sin(2*pi*tt(1:Nt)*freq);
   tone(1:NsamRamp) = tone(1:NsamRamp).*ramp;
   tone(end+1-NsamRamp:end) = tone(end+1-NsamRamp:end).*fliplr(ramp);
   Tone{ii} = Amp*tone(:);
end

% upload waveforms
seqplayupload(Tone,Tone);
% specify play list
SPinfo = seqplaylist(randperm(12), ones(1,12) ,randperm(12), ones(1,12));
%SPinfo = Seqplaylist([], [] ,randperm(12), ones(1,12));

%  monaural version
% SPinfo = Seqplayupload(Tone,{});
% LL = Seqplaylist(randperm(12), ones(1,12));
%LL = Seqplaylist(1:12, ones(1,12),1:12, ones(1,12));
%pause(0.2);
seqplaygo;





