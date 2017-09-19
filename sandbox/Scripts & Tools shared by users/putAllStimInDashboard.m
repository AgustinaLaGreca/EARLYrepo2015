%all available stim
allStim = stimCollection;
stimButtons = stimButtonTiling;
toAdd = {[]};
nrToAdd = 1;

for n = 1:size(allStim, 2)
    
    if isempty(findstr(allStim{n}, cell2words(stimButtons)))
       toAdd{ nrToAdd } = [allStim{n}];
       nrToAdd = nrToAdd +1;   
    end
    
end

       

for n = 1:nrToAdd-1
    
    for i = 1:size(stimButtons, 2)
        nrBot = size(strfind(stimButtons{i}, '/'), 2);
        if nrBot < 5
            disp(['adding: ' toAdd{n}]); 
            stimButtons{i} = [stimButtons{i}, '/', toAdd{n}];
            break;
        else
            if i == size(stimButtons, 2) && n < nrToAdd
                stimButtons{i+1} = [ toAdd{n} ];
                disp(['adding: ' toAdd{n}]);
                break;
            end
        end
    end
    
end
    
    
%stimButtonTiling(stimButtons);

    