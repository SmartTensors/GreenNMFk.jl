function [sol,normF,lb,ub, AA] = calculations_nmf(number_of_sources,nd,Nsim,aa,xD,t0,time,S,numT)                                                 
%% Calculation of the soulution with j sources 

    %j     = number_of_sources; % This loop is on the number of the unknown sources kk
    sol   = zeros(Nsim,3*number_of_sources+3);% Each solution(original source) is with a structure: [Ai Xi Yi Ux Dx Dy] i=1,2,...kk;
    normF = zeros(Nsim,1);% norm F is the OLSQ difference between the Observatio and the reconstruction
    
    
  % Defining the function that we will minimize Funf = Sum_i(MixFn(i)-Sum_k(Sources(i,k)))^2      
       for i=1:nd
           
          if  number_of_sources == 1
          
              Mixfn = @(x) source(time, x(4:6), xD(i,:), x(1:2), t0, x(3));  
          
          else
                for d=1:number_of_sources
                    if d == 1 
                        Mixfn   = @(x) source(time, x(4:6), xD(i,:), x(1:2), t0, x(3)); 
                    else
                        mixfun2 = @(x) source(time, x(d*3+1:d*3+3), xD(i,:), x(1:2), t0, x(3));
                        Mixfn   = @(x) Mixfn(x) + mixfun2(x);
                    end
                end
          end    
                 
            if i==1  

                funF = @(x) ([Mixfn(x)  zeros(1, (nd-1)*numT)]- S(i,:));
            else
                fun2 = @(x) ([ zeros(1, (i-1)*numT) Mixfn(x)  zeros(1, (nd-i)*numT)]- S(i,:));

                funF = @(x) funF(x) + fun2(x);
            end


        end
      
     % Defining the lower and upper boundary for the minimization
       lb = [0   0     0];  % lower boundary [Ux Dx Dy]
       ub = [1   1     1];  % upper boundary [Ux Dx Dy]
     % This loop is on the number of the sources we investigate 
     % we need limits for all sources (ampl and coord)
      for jj = 1:number_of_sources;
          lb = [lb 0  -aa -aa];% General lower boundary [Ux Dx Dy A X Y]
          ub = [ub 1.5 aa  aa];% General upper boundary [Ux Dx Dy A X Y]
      end
    % The norm of the observsational matrix/vector
       AA = 0;SS = 0;for i = 1:nd;SS = S(i,:).^2;AA = AA+sum(SS);end
     
        options      = optimset('MaxFunEvals', 3000);
        initCON      = zeros(Nsim,3*number_of_sources+3);
 %% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx       
         
         for k=1:Nsim % This loop iterats the NMF runs
         x_init =[rand(1) 2*aa*(0.5 - rand(1, 2))];% the IC random [A X Y]
           for d = 1:number_of_sources
         x_init = [x_init  rand() 2*aa*(0.5 - rand(1, 2))]; % the size is 3*number_of_sources+3 for the IC  
           end
         initCON(k,:) = x_init;
         end
       
        
  
    
    parfor k = 1:Nsim % This loop is iterating the NMF runs
         [sol(k,:), normF(k)] = lsqnonlin(funF,initCON(k,:),lb,ub,options);
    end
    
     
     
        normF = sqrt(normF./AA).*100;
    
            
    file_name1 = sprintf('./Results/Results_%ddet_%dsources.mat',nd, number_of_sources);
    save(file_name1, 'sol', 'normF', 'S', 'lb', 'ub', 'AA'); 
     
    
    
 end

       