function JLstimlist(Gerbil);
% JLstimlist - display stimulus list of JL exp
%    JLstimlist(iGerbil) opens teh stimulus list of the experiment of
%    gerbil # iGerbil.
%
%    See also gerbilName.

if isnumeric(Gerbil),
    Gerbil = gerbilName(Gerbil);
elseif ~isequal('RG', upper(Gerbil(1:2))),
    Gerbil = ['RG' Gerbil];
end
FN =fullfile(TKdatadir('JL'), 'StimulusLists', [ Gerbil '_SeqList.txt']);
winopen(FN);


