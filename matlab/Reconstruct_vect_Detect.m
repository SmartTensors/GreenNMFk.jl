function [Sf, Comp, Dr, Det] = CompRes(Cent,xD,Solution,t0, numT,noise, time, S )
close all
As = Cent(:,1);
Xs = Cent;
D = Solution(1,end-1:end);
u = Solution(1,end-2);

[Sf,XFf]            = initial_conditions(As,Xs,xD,D,t0,u,numT,noise,time);

Comp = sum((Sf.^2 - S.^2),2);



for i=1:size(xD,1)
    a = (i-1)*80+1;
    b = 80*i;
    Det(i,:) = S(i,a:b);
    Dr(i,:)  = Sf(i,a:b);
end



figure
g=1:2:2*size(xD,1);
for i = 1:size(xD,1)
   subplot(size(xD,1),1,i)
   plot(Det(i,:),'or');
   hold;
   %subplot(2*size(xD,2),2,g(i)+1)
   plot(Dr(i,:),'k');
end

