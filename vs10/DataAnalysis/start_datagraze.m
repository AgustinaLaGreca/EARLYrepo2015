function [] = start_datagraze()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    global datagraze
    
    if datagraze.disabled == 0
        try
            if datagraze.active == 0
                proc = System.Diagnostics.Process();
                proc.StartInfo.FileName = 'matlab';
                exp_name = current(experiment);
                exp_name = exp_name.ID.Name;
                proc.StartInfo.Arguments =  [' -nosplash ' ...
                    '-sd c:\EARLYrepo2015\vs10\dataAnalysis -r datagraze(''' exp_name ''')'];
                proc.Start();
                datagraze.active = 1;
            end
        catch
            warning(['Couldn''t start datatraze, start it manually by starting ' ...
                'a new Early matlab session and running the function ' ...
                'datagraze(ExpName) where ExpName is a String variable' ...
                ' containing the Experiment name. For expample datagraze(''a0123'')']);
            datagraze.active = 1;
        end
    end
end

