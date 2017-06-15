function [Solution, VectIndex, Cent, reconstr, mean_savg, number_of_clust_sim] = clustering_the_solutionsBSAF(number_of_sources,nd,sol,normF)
                                                                                

%K     = number_of_sources;
normF = sqrt(normF./100)*100;
    
    minsil_old = -2;
    steps      = 10;
    quants     = quantile(normF, linspace(0.2, 0.01, steps)); %start with half the data
%% Calculating the silhuettes in the different quantiles 
  u_index   = find(quants<0.1)
  u_element = max(quants(u_index));
  u = find(quants ==u_element);
  ind  = find(normF < quants(u));
  sol1 = sol(ind,:);
  
if size(sol1,1) => 5 
for u= 10:steps
       
       ind  = find(normF < quants(u));
       sol1 = sol(ind,:);
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    %            1  2  3  4  5  6   7  8  9  10 11 12
    %Each sol for three sources = [Ux Dx Dy A1 X1 Y1 A2 X2 Y2 A3 X3 Y3]
    %Each sol for two sources   = [Ux Dx Dy A1 X1 Y1 A2 X2 Y2]
    %Each sol for one source    = [Ux Dx Dy A1 X1 Y1]
       justSource = sol1(:,4:end);
       Sources3D  = zeros( 3, size(justSource,2)/3, size(justSource, 1) );

    for J = 1:size(sol1,1)
        hold = justSource(J,:);

            for kJ = 1:size(justSource,2)/3


                if kJ == 1
                    col = hold(kJ*3-2:kJ*3);
                else
                    col =[col; hold(kJ*3-2:kJ*3)];

                end

            end
        
        Sources3D(:,:,J) = col';

        if J == 1
            colSources = col;
        else  
            colSources = [colSources ; col];
        end

    end
    




    [~, VectIndex, Cent ] =  LudmilCluster( Sources3D, colSources, 100);
     
   
    [ss,h] = silhouette(colSources,VectIndex);
    savg   = grpstats(ss,VectIndex);
    minsil = savg;%min(savg);
    
    
end

number_of_clust_sim = u;

idx1     = VectIndex==1;
avg_sol  = mean(sol1(VectIndex(idx1),:));
Solution = zeros(number_of_sources,6);

    for i = 1:number_of_sources
    Solution(i,:)= [Cent(i,:) avg_sol(1:3)];
    end
    
reconstr  = mean(normF( ind));
mean_savg = mean(savg);

file_name1 = sprintf('./Results/Solution_%ddet_%dsources.mat',nd, number_of_sources);
save(file_name1, 'Solution', 'VectIndex', 'Cent', 'ss', 'savg', 'reconstr', 'mean_savg','number_of_clust_sim');   

end
 