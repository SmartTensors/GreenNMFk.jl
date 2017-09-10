
% -- Set up variables -------------------------
number_of_sources = 1;
Nsim = 100;
aa = 1;
t0 = -10;
numT = 80;	
As = [0.5; 0.5; 0.5; 0.5];
D = [0.005 0.00125];
u = 0.05;
numT = 80;
noise = 0E-3;
time = linspace(0, 20, numT);
xD = [0.0 0.0; -0.5 -0.5; -0.5 0.5; 0.5 0.5;0.5 -0.5;0.0 0.5;0.0 -0.5;-0.5 0.0; 0.5 0.0];
nd = length(xD);

Xs = ones(length(As),3);
Xn = [[-0.3; -0.4] [0.4; -0.3] [-0.1; 0.25] [-0.3; 0.65]];
for k = 1:size(Xs,1)
	Xs(k,:) = [As(k) Xn(1,k) Xn(2,k)];
end


% -- Run core functions -----------------------
tic;
[S, XF] = initial_conditions(As,Xs,xD,D,t0,u,numT,noise,time);
tc1 = toc;
[sol,normF,lb,ub, AA,sol_real, normF_real, normF1, sol_all, normF_abs, Qyes] = calculations_nmf_v02(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT);
tc2 = toc - tc1;
%[solution, vect_index, cent, reconstr, mean_savg, number_of_clust_sim] = clustering_the_solutions(number_of_sources+1, nd, sol_real, normF_real, Qyes);
tc3 = toc - tc2;

% -- Save results -----------------------------
save('./Results/initial_conditions.mat','S','XF')
save('./Results/calculations_nmf.mat','sol','normF','lb','ub','AA','sol_real','normF_real','normF1','sol_all','normF_abs','Qyes')
%save('./Results/clustering.mat','solution','vect_index','cent','reconstr','mean_savg','number_of_clust_sim')

disp("Time to compute initial conditions: " + tc1)
disp("Time to compute calculations_nmf:   " + tc2)
disp("Time to compute clustering:         " + tc3)
