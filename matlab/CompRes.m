function [Sf, Comp, Dr, Det, Wf] = CompRes(Cent,xD,Solution,t0, numT,noise, S )
%close all
time = linspace(0, 20, numT);
if  (size(Solution,1) == 1)    
  As =Solution(4);
  Xs = Solution(4:6);
  D  = Solution(1:2);
  u  = Solution(3);
else
 As = Cent(:,1);
 Xs = Cent;
 D = Solution(1,end-2:end-1);
 u = Solution(1,end);
end

[Sf,XFf, Wf]            = initial_conditions_2(As,Xs,xD,D,t0,u,numT,noise,time);

figure(2)
bar(Wf,'stack')

Comp = sum((Sf.^2 - S.^2),2);



for i=1:size(xD,1)
    a = (i-1)*80+1;
    b = 80*i;
    Det(i,:) = S(i,a:b);
    Dr(i,:)  = Sf(i,a:b);
end



figure(3)
g=1:2:2*size(xD,1);
for i = 1:size(xD,1)
   subplot(size(xD,1),1,i)
   plot(Det(i,:),'or');
   hold;
   %subplot(2*size(xD,2),2,g(i)+1)
   plot(Dr(i,:),'k');
end

