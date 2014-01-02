function y=simulationPDE
NCLOCK=fix(clock)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS OF THE PDE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbday=400;    % in days                           % number of days for the simulation

amax=90;        % in days 
time_step=1/3;  % in days
h=time_step;
nb_age_step=round(amax/time_step);

nb_time_step=round(nbday/time_step);     % number of time step per shift


% Parameters concerning the hospital 

nbh=100;                % number of healthcare workers 
nbp=400;                % number of patients


AU=5;                 % Average Length of Stay (LOS) for patients 'u'
AN=14;                % Average Length of Stay (LOS) for patients 'n'
AR=28;                % Average Length of Stay (LOS) for patients 'r'

% Parameters concerning infection, contamination, visits 

AV=60;                  % Average time of Visit in minutes
AV=AV/(60*24); % in day unit   

AC=60;                  % Average time of contamination of hcws in minutes  
AC=AC/(60*24); % in day unit  

PC=0.4;                 % Probability of contamination during 1 visit
PI=0.06;                % Probability of infection during 1 visit
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial patients distribution 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
PNN0=7;                 % precentage of patients 'n' at time t=0
PRS0=0;                 % precentage of patients 'rs' at time t=0
PRR0=0;                 % precentage of patients 'rr' at time t=0
PNR0=3;                 % precentage of patients 'nr' at time t=0
                   
PNN0=PNN0/100;
PRS0=PRS0/100;
PRR0=PRR0/100;
PNR0=PNR0/100;
PUU0=1-(PNN0+PRS0+PRR0+PNR0);

AI=7;               % Average infection age at time t=0 in days

NUI=1/AI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters of the bacterial load development for 1 patient
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bt=3;       %beginning of treatment 
et=21;       %end of treatment 

TN=10^11;    % Infectiousness threshold with respect to N-strain     
TR=10^11;    % Infectiousness threshold with respect to R-strain 
TREC=10^5;   % threshold of recovery 

tau=0.001;    % recombinasion rate 
kappa=10^15;  %
gamma=10^-25; % Reversion rate 

% Growth rates of bacteria when there is no treatment 
betaN1=4;
betaR1=2;
% Growth rates of bacteria during the treatment 
betaN2=-3;
betaR2=2;
% Growth rate of the bacterial load is <100. 
betaN3=-3;
betaR3=-3;


% END OF THE PARAMETERS LIST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compute the parameter of the PDE model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NUV=1/AV;

NUC=1/AC;

ALPHA=PC*NUV;

NUN=1/AN;

NUR=1/AR;

BETAV=nbh/nbp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compute the infectiousness functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for patient 'n'  
patient=newpatientN;
V=zeros(1,2);
V=[10^7 0]; 
gn_n=zeros(1,nb_age_step+1);
for i=1:nb_age_step+1
    t(i)=h*(i-1);
    if (V(1)<TN)
        gn_n(i)=0;
    else 
        gn_n(i)=1;
    end
        
    patient=evolutiondosepatient(h,bt,et,TN,TR,TREC,tau,kappa,gamma,betaN1,betaR1,betaN2,betaR2,betaN3,betaR3,patient);
    V(1)=patient{3}(1);
    V(2)=patient{3}(2);
    
end

% for patient 'rs'  
patient=newpatientRS;
V=zeros(1,2);
V=[10^10 5*10^6];
grs_n=zeros(1,nb_age_step+1);
grs_nr=zeros(1,nb_age_step+1);
grs_r=zeros(1,nb_age_step+1);
for i=1:nb_age_step+1
    if (V(1)<TN & V(2)<TR)
        grs_n(i)=0;
        grs_nr(i)=0;
        grs_r(i)=0;
    else 
        if (V(1)>=TN & V(2)>=TR)
            grs_n(i)=0;
            grs_nr(i)=1;
            grs_r(i)=0;
        else 
            if (V(1)>=TN & V(2)<TR)
                grs_n(i)=1;
                grs_nr(i)=0;
                grs_r(i)=0;
            else
                grs_n(i)=0;
                grs_nr(i)=0;
                grs_r(i)=1;
            end
        end
    end
    patient=evolutiondosepatient(h,bt,et,TN,TR,TREC,tau,kappa,gamma,betaN1,betaR1,betaN2,betaR2,betaN3,betaR3,patient);
    V(1)=patient{3}(1);
    V(2)=patient{3}(2);
    
end

% for patient 'rr'  
patient=newpatientRR;
V=zeros(1,2);
V=[0 5*10^6];
grr_n=zeros(1,nb_age_step+1);
grr_nr=zeros(1,nb_age_step+1);
grr_r=zeros(1,nb_age_step+1);
for i=1:nb_age_step+1
    if (V(1)<TN & V(2)<TR)
        grr_n(i)=0;
        grr_nr(i)=0;
        grr_r(i)=0;
    else 
        if (V(1)>=TN & V(2)>=TR)
            grr_n(i)=0;
            grr_nr(i)=1;
            grr_r(i)=0;
        else 
            if (V(1)>=TN & V(2)<TR)
                grr_n(i)=1;
                grr_nr(i)=0;
                grr_r(i)=0;
            else
                grr_n(i)=0;
                grr_nr(i)=0;
                grr_r(i)=1;
            end
        end
    end
    patient=evolutiondosepatient(h,bt,et,TN,TR,TREC,tau,kappa,gamma,betaN1,betaR1,betaN2,betaR2,betaN3,betaR3,patient);
    V(1)=patient{3}(1);
    V(2)=patient{3}(2);
    
end




% for patient 'nr'  
patient=newpatientR;
V=zeros(1,2);
V=[10^7 5*10^6]; 
gnr_n=zeros(1,nb_age_step+1);
gnr_nr=zeros(1,nb_age_step+1);
gnr_r=zeros(1,nb_age_step+1);
for i=1:nb_age_step+1
    
    if (V(1)<TN & V(2)<TR)
        gnr_n(i)=0;
        gnr_nr(i)=0;
        gnr_r(i)=0;
    else 
        if (V(1)>=TN & V(2)>=TR)
            gnr_n(i)=0;
            gnr_nr(i)=1;
            gnr_r(i)=0;
        else 
            if (V(1)>=TN & V(2)<TR)
                gnr_n(i)=1;
                gnr_nr(i)=0;
                gnr_r(i)=0;
            else
                gnr_n(i)=0;
                gnr_nr(i)=0;
                gnr_r(i)=1;
            end
        end
    end
    patient=evolutiondosepatient(h,bt,et,TN,TR,TREC,tau,kappa,gamma,betaN1,betaR1,betaN2,betaR2,betaN3,betaR3,patient);
    V(1)=patient{3}(1);
    V(2)=patient{3}(2);
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compute the initial distributions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PUU(1)=PUU0;

for i=1:nb_age_step+1
DPNN(1,i)=PNN0*NUI*exp(-NUI*(i-1)*h);
DPRS(1,i)=PRS0*NUI*exp(-NUI*(i-1)*h);
DPRR(1,i)=PRR0*NUI*exp(-NUI*(i-1)*h);
DPNR(1,i)=PNR0*NUI*exp(-NUI*(i-1)*h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=zeros(1,nb_time_step+1);
PN=zeros(1,nb_time_step+1);
PRS=zeros(1,nb_time_step+1);
PRR=zeros(1,nb_time_step+1);
PNR=zeros(1,nb_time_step+1);
for j=1:nb_time_step
    time(j+1)=h*j;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We compute the probabilities to be infectious PIN PINR PIR
    % and the probabilities PN(t) and PRR(t),PRS(t), PNR(t) to belong to 'n' and 'r' 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for i=1:nb_age_step+1
        DN=gn_n.*DPNN(j,:)+grs_n.*DPRS(j,:)+grr_n.*DPRR(j,:)+gnr_n.*DPNR(j,:);
        DNR=grs_nr.*DPRS(j,:)+grr_nr.*DPRR(j,:)+gnr_nr.*DPNR(j,:);
        DR=grs_r.*DPRS(j,:)+grr_r.*DPRR(j,:)+gnr_r.*DPNR(j,:);
        %end
    n=nb_age_step/2;
    SN1=0; 
    SN2=0;
    SNR1=0;
    SNR2=0;
    SR1=0;
    SR2=0;
    WN1=0;
    WN2=0;
    WRS1=0;
    WRS2=0;
    WRR1=0;
    WRR2=0;
    WNR1=0;
    WNR2=0;
    for i=1:n-1
        SN1=SN1+DN(2*i+1);
        SNR1=SNR1+DNR(2*i+1);
        SR1=SR1+DR(2*i+1);
        WN1=WN1+DPNN(j,2*i+1);
        WRS1=WRS1+DPRS(j,2*i+1);
        WRR1=WRR1+DPRR(j,2*i+1);
        WNR1=WNR1+DPNR(j,2*i+1);
    end
    for i=1:n
        SN2=SN2+DN(2*i);
        SNR2=SNR2+DNR(2*i);
        SR2=SR2+DR(2*i);
        WN2=WN2+DPNN(j,2*i);
        WRS2=WRS2+DPRS(j,2*i);
        WRR2=WRR2+DPRR(j,2*i);
        WNR2=WNR2+DPNR(j,2*i);
    end
    % probabilities to be infectious with respect to n nr and r
    PIN=(h/3)*(DN(1)+2*SN1+4*SN2+DN(nb_age_step+1));
    PINR=(h/3)*(DNR(1)+2*SNR1+4*SNR2+DNR(nb_age_step+1));
    PIR=(h/3)*(DR(1)+2*SR1+4*SR2+DR(nb_age_step+1));
    
    
    PN(j)=(h/3)*(DPNN(j,1)+2*WN1+4*WN2+DPNN(j,nb_age_step+1));
    PRS(j)=(h/3)*(DPRS(j,1)+2*WRS1+4*WRS2+DPRS(j,nb_age_step+1));
    PRR(j)=(h/3)*(DPRR(j,1)+2*WRR1+4*WRR2+DPRR(j,nb_age_step+1));
    PNR(j)=(h/3)*(DPNR(j,1)+2*WNR1+4*WNR2+DPNR(j,nb_age_step+1));
    PU(j)=1-(PN(j)+PRS(j)+PRR(j)+PNR(j));
    
    HU=NUC/(NUC+ALPHA*(PIN+PINR+PIR));
    HN=(ALPHA*PIN/(NUC+ALPHA*(PIR+PINR)))*HU;
    HR=(ALPHA*PIR/(NUC+ALPHA*(PIN+PINR)))*HU;
    HNR=ALPHA*((PINR+PIR)*HN+PINR*HU+(PIN+PINR)*HR)/NUC;
        
    t0=h*(j-1);
    t2=h*j;
    t1=(t0+t2)/2;
    T=[t0 t1 t2];
    patient0=zeros(1,5);
    patient0=[PU(j) PN(j) PRS(j) PRR(j) PNR(j)];
    
    [t,x]=ode45(@EDOpatientinfection,T,patient0,[],NUN,NUR,NUV,BETAV,PI,HN,HNR,HR);
    PU(j+1)=x(3,1);
    PN(j+1)=x(3,2);
    PRS(j+1)=x(3,3);
    PRR(j+1)=x(3,4);
    PNR(j+1)=x(3,5);
     
    
    
    AUX2=(NUV*BETAV*(HNR+HR)+NUN);
    for i=nb_age_step+1:-1:2
      DPNN(j+1,i)=exp(-h*AUX2)*DPNN(j,i-1);
      DPRS(j+1,i)=exp(-h*NUR)*DPRS(j,i-1);
      DPRR(j+1,i)=exp(-h*NUR)*DPRR(j,i-1);
      DPNR(j+1,i)=exp(-h*NUR)*DPNR(j,i-1);
    end
    
    DPNN(j+1,1)=NUV*BETAV*HN*PI*PU(j+1);
    DPRS(j+1,1)=NUV*BETAV*(HNR+HR)*PI*PN(j+1);
    DPRR(j+1,1)=NUV*BETAV*HR*PI*PU(j+1);
    DPNR(j+1,1)=NUV*BETAV*HNR*PI*PU(j+1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compute the R0N and R0R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEXPN=zeros(1,nb_age_step+1);
DEXPR=zeros(1,nb_age_step+1);
for i=1:nb_age_step+1
DEXPN(i)=exp(-NUN*(i-1)*h);
DEXPR(i)=exp(-NUR*(i-1)*h);
end
DN=gn_n.*DEXPN;
DR_R=(grr_r).*DEXPR;
DNR_R=(gnr_r).*DEXPR;
DR_NR=(grr_nr).*DEXPR;
DNR_NR=(gnr_nr).*DEXPR;

SN1=0; 
SN2=0;
SR_R1=0;
SR_R2=0;
SNR_R1=0;
SNR_R2=0;
SR_NR1=0;
SR_NR2=0;
SNR_NR1=0;
SNR_NR2=0;

for i=1:n-1
    SN1=SN1+DN(2*i+1);
    
    SR_R1=SR_R1+DR_R(2*i+1);
    SNR_R1=SNR_R1+DNR_R(2*i+1);
    
    SR_NR1=SR_NR1+DR_NR(2*i+1);
    SNR_NR1=SNR_NR1+DNR_NR(2*i+1);
    
   
end
for i=1:n
    SN2=SN2+DN(2*i);
  %  SRR2=SRR2+DRR(2*i);
  %  SNR2=SNR2+DNR(2*i);
    
    SR_R2=SR_R1+DR_R(2*i);
    SNR_R2=SNR_R1+DNR_R(2*i);
    
    SR_NR2=SR_NR2+DR_NR(2*i);
    SNR_NR2=SNR_NR2+DNR_NR(2*i);
end
% probabilities to be infectious with respect to n nr and r
SN=(h/3)*(DN(1)+2*SN1+4*SN2+DN(nb_age_step+1));

SR_R=(h/3)*(DR_R(1)+2*SR_R1+4*SR_R2+DR_R(nb_age_step+1));
SNR_R=(h/3)*(DNR_R(1)+2*SNR_R1+4*SNR_R2+DNR_R(nb_age_step+1));

SR_NR=(h/3)*(DR_NR(1)+2*SR_NR1+4*SR_NR2+DR_NR(nb_age_step+1));

SNR_NR=(h/3)*(DNR_NR(1)+2*SNR_NR1+4*SNR_NR2+DNR_NR(nb_age_step+1));



R0N=NUV*BETAV*PI*ALPHA*SN/NUC; 

MAT=[SR_R SNR_R; SR_NR SNR_NR];

aux=eig(MAT);
R0R=NUV*BETAV*PI*ALPHA*max(aux)/NUC;



%R0RR=NUV*BETAV*PI*ALPHA*SRR/NUC; 
%R0NR=NUV*BETAV*PI*ALPHA*SNR/NUC; 








PN=100*PN;
PRS=100*PRS;
PRR=100*PRR;
PNR=100*PNR;
PR=PRS+PRR+PNR;
PT=PN+PR;
plot(time,PT,time,PN,time,PR,time,PNR,time,PRS,time,PRR,'LineWidth',1.5);
h=legend('% of N and R','% of N','% of R','% of NR','% of RS','% of RR');
xlabel('days');
ylabel('Percentage of patients');
AV=AV*(60*24);
AC=AC*(60*24);
PNN0=PNN0*100;
PRS0=PRS0*100;
PRR0=PRR0*100;
PNR0=PNR0*100;

title([' AV=',num2str(AV),' mns, AC=',num2str(AC),' mns, bt=',num2str(bt),', et=',num2str(et),', R0N=',num2str(R0N,3),', R0R=',num2str(R0R,3)]);
axis([-Inf Inf 0 100])
NCLOCK=fix(clock)