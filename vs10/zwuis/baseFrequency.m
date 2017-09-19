function [Fbase, NsamBuf, BufDur] = baseFrequency(Fsam, ApproxBufDur)
% baseFrequency - evaluate "frequency quantum" compatible with sample rate
%   [Fbase, NsamBuf, BufDur] = baseFrequency(Fsam, ApproxBufDur) evaluates
%   the lowest nonzero frequency Fbase that exactly fits in a buffer with
%   approximate duration ApproxBufDur. All other fitting frequencies are
%   multiples of Fbase.
%
%   Inputs
%          Fsam: sample frequency in kHz
%  ApproxBufDur: approximate circular buffer duration in ms
%
%   Outputs
%         Fbase: base frequency in Hz (!!) as described above
%       NsamBuf: # samples in circular buffer
%        BufDur: exact duration of circular buffer in ms

SamPeriod = 1./Fsam; % sample period in ms
NsamBuf = round(ApproxBufDur./SamPeriod); % nearest multiple of SamPeriod
BufDur = NsamBuf.*SamPeriod; % exact buffer dur in ms
Fbase = 1e3./BufDur; % freq in Hz corresponding to BufDur




