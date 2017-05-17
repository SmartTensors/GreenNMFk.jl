function createFigure2(Wf,nd,nopt)
%CREATEFIGURE(YMATRIX1)
%  YMATRIX1:  bar matrix data

%  Auto-generated by MATLAB on 23-Sep-2016 12:57:54

% Create figure

sources = cell(nopt, 1);
for i=1:nd
   sources{i} = strcat('S', num2str(i));
end


detectors = cell(nd, 1);
for i=1:nd
   detectors{i} = strcat('D', num2str(i));
end

TICKs = 1:1:nd;

figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',detectors,...
    'XTick',TICKs,...
    'FontSize',16);
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[-0.4 10]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to bar
bar1 = bar(Wf,'BarLayout','stacked','Parent',axes1);
for i=1:nopt
set(bar1(i),'DisplayName',sources{i});
end
% Create xlabel
xlabel('Sensors','FontWeight','bold','FontSize',30);

% Create ylabel
ylabel({'Total contribution of the sources in the sensor mixtures [%]'},...
    'FontWeight','bold',...
    'FontSize',24);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1],...
    'Position',[0.0318968739637212 0.435648589732057 0.0496158770806658 0.110636645962733]);

