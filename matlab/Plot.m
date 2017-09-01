clear
Path = './Results5d3s_Good/';
nd = str2num(Path(10));
ns = str2num(Path(12));

SILLAVG = zeros(length(nd),1);
RECON = zeros(length(nd),1);

for i = 1:nd
    file = sprintf('Solution_%ddet_%dsources.mat',nd, i);
    file = [Path file];
    
    load (file);
    SILLAVG(i,1) = mean_savg;
    
    if i ==1
        RECON(i,1) = reconstr1/37.0271;
        SILLAVG(i,1) = 1;
    else
         RECON(i,1) = reconstr;
    end
    
end



[hl1, hl2, x, y1, y2] = PlotResultsBSA(RECON, SILLAVG, nd)

