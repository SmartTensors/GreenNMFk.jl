clear all
close all

max_number_of_sources = 9;
Nsim                  = 100;

% The # of the time points
numT = 80;
time = linspace(0, 20, numT);

u     = 0.05;% [km/year]; %flow speed
D     = [0.005  0.00125] ; % difusion coefficient [km^2/year; Dx=Dispersivity_long*u;Dx=Dispersivity_transv*u;
t0    = -10; % the initial time of the sources
noise = 0*10^(-3); %the noise strenght

% The amplitudes of the real sources
As    = [0.5; 0.5; 0.5; 0.5]; % Write desired amplitudes in a vector [A1 A2...]

% The real number of sources
ns    = length(As);

% Xn -> the real sources's coordinates
Xn = zeros(2,length(As));
aa = 1;% multipl for the initial random conditions

%The coordinates of the sources
Xn = [[-0.3;   -0.4] [0.4; -0.3] [ -0.1; 0.25] [ -0.3; 0.65]];

Xs = ones(length(As), 3);

% Ordering the matrix of the sources: [A X Y]
for k = 1:size(Xs,1)
 Xs(k,:) = [As(k) Xn(1,k) Xn(2,k)];
end

xd1 =  [0.0  0.0]; % position of the first detector
xd2 =  [-0.5 -0.5] ; % position of the second detector
xd3 =  [-0.5  0.5]; % position of the third detector
xd4 =  [0.5   0.5]; % position of the fourth detector
xd5 =  [0.5  -0.5];
xd6 =  [0.0   0.5] ; % position of the second detector
xd7 =  [0.0  -0.5]; % position of the third detector
xd8 =  [-0.5  0.0]; % position of the fourth detector
xd9 =  [0.5   0.0];

xD = [xd1;xd2;xd3;xd4; xd5;xd6;xd7;xd8; xd9];
nd = length(xD);

aa = 1;% the length of the interaval for random IC

clear normF sol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S,XF]            = initial_conditions(As,Xs,xD,D,t0,u,numT,noise,time);
fprintf('\nSaving init. conditions: xtrue_%ddet_%dsources.mat\n',nd, length(As))

number_of_sources = 1;
Nsim = 5

calculations_nmf_v02(number_of_sources,nd,Nsim,aa,xD,t0,time,S,numT);
fprintf('\nSaving results: Results_%ddet_%dsources.mat\n', nd, number_of_sources)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

