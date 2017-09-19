function [ Pres ] = addStutterToPres( Pres )
%addStutterToPres Make the presentation field stutter compatible

Npres = Pres.Npres;
Ncond = Pres.Ncond;

Pres.iCond = [Pres.iCond(1); 1; Pres.iCond(2:end)];
Pres.iRep = [Pres.iRep(1); 1; Pres.iRep(2:end)];
Pres.NsamPres = [Pres.NsamPres(1); Pres.NsamPres(2); Pres.NsamPres(2:end)]; 
Pres.SamOffset = [Pres.SamOffset(1:2);Pres.SamOffset(3); Pres.NsamPres(2)+Pres.SamOffset(3:end)];
Pres.PresOnset = 1e3*Pres.SamOffset/Pres.Fsam;
Pres.PresDur = [Pres.PresDur(1); Pres.PresDur(2); Pres.PresDur(2:end)];

% Update the X and or Y field of the Presentation Struct
if ~isempty(Pres.X)
    Pres.X.PlotVal = [Pres.X.PlotVal];
    Pres.X.Ncond = Pres.X.Ncond;
end

if ~isempty(Pres.Y)
    Pres.Y.PlotVal = [Pres.Y.PlotVal];
    Pres.Y.Ncond = Pres.Y.Ncond;
end

Pres.Ncond = Pres.Ncond;
Pres.Npres = Pres.Npres+1;

Pres.TotDur = Pres.TotDur + Pres.PresDur(2);
Pres.TotNsamPlay = Pres.TotNsamPlay+Pres.NsamPres(2);

end

