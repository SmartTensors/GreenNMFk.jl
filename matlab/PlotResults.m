function [hl1, hl2] = PlotResults(RECON, SILLAVG, number_of_sources)
K = number_of_sources;

figure
x1 = [1:1:K];
y1 = RECON;% 4.*cos(x1)./(x1+2);
x2 = [1:1:K];
y2 = SILLAVG;%x2.^2./x2.^3;

%Using low-level line and axes routines allows you to superimpose objects easily. Plot the first data, making the color of the line and the corresponding x- and y-axis the same to more easily associate them.
hl1 = line(x1,y1,'Color','r', 'Marker', '.', 'MarkerSize', 30);
ax1 = gca;
set(ax1,'XColor','r','YColor','r')



%Next, create another axes at the same location as the first, placing the x-axis on top and the y-axis on the right. Set the axes Color to none to allow the first axes to be visible and color code the x- and y-axis to match the data.

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k', 'title', 'Yes');
title('Results');
xlabel(ax1, 'The number of original sources');
ylabel(ax1 ,'Reconstruction Error [a.u.]');
ylabel(ax2 , 'Silhouettes mean vale');


%Draw the second set of data in the same color as the x- and y-axis.
hl2 = line(x2,y2,'Color','k','Parent',ax2, 'Marker', '.', 'MarkerSize', 30);

end