tic;

% delete(gcp)
% parpool('local')
% myCluster = parcluster('local');

matlabpool open


%% Initial conditions
clear all
close all

max_number_of_sources = 9;
Nsim                  = 100;
% The # of the time points
numT = 80; 
time = linspace(0, 20, numT);

% The flow velocity (it is [0.05 0] in km/year



u     = 0.05;% [km/year]; %flow speed
% Dispersivity_long = 0.1;%[km]
% Dispersivity_long = 0.025;%[km]
D     = [0.005  0.00125] ; % difusion coefficient [km^2/year; Dx=Dispersivity_long*u;Dx=Dispersivity_transv*u;
t0    = -10; % the initial time of the sources
noise = 1*10^(-4); %the noise strenght
% The amplitudes of the real sources
%As    = [0.5]; %; 0.7];% 0.3]; % Write desired amplitudes in a vector [A1 A2...]
As = 0.3;
% The real number of sources
ns    = length(As);

% Xn -> the real sources's coordinates
Xn = zeros(2,length(As));
aa = 1;% multipl for the initial random conditions
%The coordinates of the sources
%Xn = [-0.9;   -0.8];% [-0.1; -0.2]];% [ -0.2; 0.6]];
Xn = [-0.1; -0.2];
Xs = ones(length(As), 3);
% Ordering the matrix of the sources: [A X Y]
for k = 1:size(Xs,1)
 Xs(k,:) = [As(k) Xn(1,k) Xn(2,k)];
end
%  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% xd -> the positions of the detectors
xd1 =  [0  0]; % position of the first detector
%xd2 =  [-0.5 -0.5] ; % position of the second detector
%xd3 =  [-0.5  0.5]; % position of the third detector
xd4 =  [0.5   0.5]; % position of the fourth detector
xd5 =  [0.5  -0.5]; 
% xd6 =  [ 0  0.5]; 
% xd7 =  [ 0 -0.5]; 
% xd8 =  [-0.5  0];
% xd9 =  [ 0.5  0]; 
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

xD = [xd1;xd4; xd5];
%xD = [xd1; xd2; xd3; xd4; xd5];
% The number of detectors
nd = length(xD);

aa = 1;% the length of the interaval for random IC

clear normF sol
%% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    [S,XF]            = initial_conditions(As,Xs,xD,D,t0,u,numT,noise,time);

    number_of_sources = 1;
    disp(number_of_sources);
    
    [sol,normF,lb,ub, AA, sol_real,normF_real, normF1, sol_all] = calculations_nmf_v02(number_of_sources,nd,Nsim,aa,xD,t0,time,S,numT);
    disp(size(sol));
    yy = quantile(normF,.25);
    reconstr1 = mean(normF(normF<yy));
    
    ind  = find(normF < yy);
    sol1 = sol(ind,:);
    
    avg_sol             = mean(sol1);
    Solution            = avg_sol;
    mean_savg           = 1;
    number_of_clust_sim = 0;

    file_name1 = sprintf('./Results/Solution_%ddet_%dsources.mat',nd, number_of_sources);
    save(file_name1, 'Solution', 'reconstr1', 'mean_savg','number_of_clust_sim');  

RECON(1)   = reconstr1;
SILL_AVG(1) = 1;


number_of_sources = 2;

for jj = 2:max_number_of_sources
  
disp(number_of_sources);

[sol,normF,lb,ub, AA,sol_real, normF_real, normF1, sol_all, normF_abs, Qyes] = calculations_nmf_v02(number_of_sources,nd,Nsim,aa,xD,t0,time,S,numT);


[Solution, VectIndex, Cent, reconstr, mean_savg, number_of_clust_sim] = clustering_the_solutions(number_of_sources,nd,sol_real,normF_real, Qyes);

RECON(number_of_sources)    = reconstr;
SILL_AVG(number_of_sources) = mean_savg;

number_of_sources = number_of_sources + 1;
end
close all
RECON = RECON/Nsim; 
file_name1 = sprintf('./Results/All_%ddet_%dsources.mat',nd, length(As));
save(file_name1, 'RECON', 'SILL_AVG'); 


x = 1:1:max_number_of_sources;
y1 = RECON;% 4.*cos(x1)./(x1+2);
y2 = SILL_AVG;

createfigureBSA(x, y1, y2)

[Sf, Comp, Dr, Det, Wf] = CompRes(Cent,xD,Solution,t0, numT,noise, time, S );
matlabpool close  
toc
