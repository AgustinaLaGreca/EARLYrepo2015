function StimStr = durationstring(D, flag);
% Dataset/durationstring - char string describing durations and reps
%   durationstring(DS) returns a string like
%     34 x 10 x 70/100 ms   
%   meaning
%     #conditions x #reps x burstdur/repdur
%   
%   If DS has two varied stimulus parameters X & Y then the strings is like
%     23 x 9 x 20 x 70/100 ms   
%   meaning
%     #X-valuess x #Y-values #reps x burstdur/repdur
%
%   durationstring(DS, 'compact') squeezes out the blanks, except the one
%   before 'ms.'
%
%   See dataset/repdur, dataset/burstdur, dataset/stimlist, dataset/disp.

[flag] = arginDefaults('flag', '');

if isequal('compact', flag),
    tstr = 'x';
else,
    tstr = ' x ';
end

Dur = burstdur(D,1);
Dur = unique(0.1*round(10*mean(Dur,1)));
DurStr = [local_str(Dur) '/' num2str(unique(0.1*round(10*mean(repdur(D),1)))) ' ms'];
Pres = D.Stim.Presentation;
if has2varparams(D),
    StimStr = [num2str(D.Stim.Ncond_XY(1)) tstr num2str(D.Stim.Ncond_XY(2)) tstr num2str(Pres.Nrep) tstr DurStr];
else,
    StimStr = [num2str(Pres.Ncond) tstr num2str(Pres.Nrep) tstr DurStr];
end
if ~isempty(D.Stopped),
    StimStr = [StimStr ' ***'];
end

%===========================
function s=local_str(x,sep);
if nargin<2, sep='|'; end
x = unique(x);
if numel(x)>1,
    s = [local_str(x(1)) sep local_str(x(2:end))];
    return;
end
if abs(x)>100,
    s = num2str(round(x));
else
    s = num2str(deciRound(x,3));
end








