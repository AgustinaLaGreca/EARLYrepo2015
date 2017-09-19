function [] = close_datagraze()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    global datagraze
    
    if datagraze.disabled == 0
        if datagraze.active == 1
            matlab_processes = System.Diagnostics.Process.GetProcessesByName('matlab');
            this_session = System.Diagnostics.Process.GetCurrentProcess();
            if matlab_processes.Length == 2 && datagraze.active==1
               mat1 = matlab_processes(1);
               mat2 = matlab_processes(2);
               if mat1.Id == this_session.Id
                   try
                       mat2.Kill()
                       datagraze.active = 0;
                   catch
                       warning('Cannot close datagraze please do it manually!');
                       datagraze.active = 0;
                   end
               elseif mat2.Id == this_session.Id
                   try
                       mat1.Kill()
                       datagraze.active = 0;
                   catch
                       warning('Cannot close datagraze please do it manually!');
                       datagraze.active = 0;
                   end
               end
               disp('Datagraze closed');
            elseif matlab_processes.Length > 2
                warning(['Cannot close datagraze because more than two matlab sessions '...
                    'are open, please close the datagraze command window!']);
                datagraze.active = 0;
            elseif matlab_processes.Length == 1
                warning('Cannot close datagraze since it isn''t open anymore!');
                datagraze.active = 0;
            end
        end
    end
end

