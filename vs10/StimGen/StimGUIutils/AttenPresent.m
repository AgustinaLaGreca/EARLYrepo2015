function ap=AttenPresent;
% AttenPresent - true if analog attenuators are present
%    AttenPresent returns True (1) if anaolog attenuators are present. 
%    Currently, this is equivalent to the presence of PA5_1 and PA5_2 in
%    the sys3devicelist, but that might change one day ;)

ap = all(ismember({'PA5_1', 'PA5_2'}, sys3devicelist));




