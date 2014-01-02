function dy=Bacterialload(t,y,tau,kappa,gamma,TREC,betaN,betaR,betaN1,betaR1)
%
% The fonction corresponds to the second member of the ODE model 
% during the treatment period
% 
% Note that we take the same values as in PNAS paper fig. 3.5 
%
% betaN=beta-
% betaR=beta+
%
%
dy=zeros(2,1);  
if (y(1)+y(2)<TREC)
    if (y(1)+y(2)==0)
        dy(1)=0;
        dy(2)=0;
    else 
        dy(1)=(-tau*y(2)/(y(1)+y(2))+betaN1-(y(1)+y(2))/kappa)*y(1)+gamma*y(2);
        dy(2)=(tau*y(1)/(y(1)+y(2))+betaR1-(y(1)+y(2))/kappa)*y(2)-gamma*y(2);
    end
else
    dy(1)=(-tau*y(2)/(y(1)+y(2))+betaN-(y(1)+y(2))/kappa)*y(1)+gamma*y(2);
    dy(2)=(tau*y(1)/(y(1)+y(2))+betaR-(y(1)+y(2))/kappa)*y(2)-gamma*y(2);
   

end 