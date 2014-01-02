function Y=simulation
NCLOCK=fix(clock)
%
% This program contains the full IBM model with evolution of the
% infectiousness of patients only at the end of each shift, but the exits
% of patients may happen at any time during the shift. 
%

%
% At the end of the program we save the figure into a file called
%  name-of-the-figure.fig
%

filename='name-of-the-figure.fig';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS OF THE IBM MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters concerning the time

nbday=1;                            % number of days for the simulation

nbs=1;                              % number of simulations (the result will be the average trajectory)

shift=8;                            % length of shift in hours

shift=shift/24;                     % length of shift in day unit

nbshift=3;                          % number of shifts per day

timestep=5;                        % time step in the shift in minutes

timestep=timestep/(24*60);          % time step in day unit 

nbtimestep=fix(shift/timestep);     % number of time step per shift

timestepinfectionage=shift;         % timestep for the age of infection
% Parameters concerning the hospital 

nbh=100;                % number of healthcare workers 
nbp=400;                % number of patients

ALSU=5;                 % Average Length of Stay (LOS) for patients 'u'
ALSN=14;                % Average Length of Stay (LOS) for patients 'n'
ALSR=28;                % Average Length of Stay (LOS) for patients 'r'

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
                  
PN=7;                    % precentage of patients 'n' at time t=0 (i.e. at the beginning of the simulation)
PRS=0;                  % precentage of patients 'rs' at time t=0
PRR=0;                  % precentage of patients 'rr' at time t=0
PNR=3;                  % precentage of patients 'nr' at time t=0
                   
PN=PN/100;
PRS=PRS/100;
PRR=PRR/100;
PNR=PNR/100;
PU=1-(PN+PRS+PRR+PNR);

AI=7;               % Average infection age at time t=0 in days

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters of the bacterial load development for 1 patient
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bt=3;        %beginning of treatment 
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

          

% END OF THE PARAMETER LIST


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          WE DEFINE RANDOM VARIABLES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXPONENTIAL RANDOM LAWS


expexitU=newexprnd(ALSU);           % exit time for patient 'u' 

expexitN=newexprnd(ALSN);           % exit time for patient 'n'    

expexitR=newexprnd(ALSR);           % exit time for patient 'r'    

expvisit=newexprnd(AV);             % length of visits

expcontamhcw=newexprnd(AC);         % Time of contamination of hcws  

expinit=newexprnd(AI);              % Initial infection-age distribution 

% UNIFORM RANDOM LAWS

random1=newrandom; % For the contamination
random2=newrandom; % For the infection 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     INITIALIZATION OF PATIENTS DISTRIBUTION AT TIME 0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%n1=round(nbp*PN);
n2=round(nbp*PRS);
n3=round(nbp*PRR);
n4=round(nbp*PNR);

n5=n1+n2;
n6=n5+n3;
n7=n6+n4;

for j=1:nbs

    for k=1:n1
        %patient{1}{j,k}=newpatientN;
        patient{1}{j,k}=cell(1,11);
        patient{1}{j,k}{1}=0;
        patient{1}{j,k}{2}='n';
        patient{1}{j,k}{3}=[10^7 0];
        patient{1}{j,k}{4}=0;
        patient{1}{j,k}{5}=0;
        patient{1}{j,k}{6}='n'; 
        patient{1}{j,k}{7}='n';
        patient{1}{j,k}{8}='ui';
        patient{1}{j,k}{9}=0;
        patient{1}{j,k}{10}=0; 
        patient{1}{j,k}{11}='f';
    end
    m=n1+1;
    for k=m:n5
        %patient{1}{j,k}=newpatientRS;
        patient{1}{j,k}=cell(1,11);
        patient{1}{j,k}{1}=0;
        patient{1}{j,k}{2}='r';
        patient{1}{j,k}{3}=[10^10 5*10^6];
        patient{1}{j,k}{4}=0;
        patient{1}{j,k}{5}=0;
        patient{1}{j,k}{6}='n'; 
        patient{1}{j,k}{7}='n';
        patient{1}{j,k}{8}='rs';
        patient{1}{j,k}{9}=0;
        patient{1}{j,k}{10}=0; 
        patient{1}{j,k}{11}='f';

    end
    m=n5+1;
    for k=m:n6
        %patient{1}{j,k}=newpatientRR;
        patient{1}{j,k}=cell(1,11);
        patient{1}{j,k}{1}=0;
        patient{1}{j,k}{2}='r';
        patient{1}{j,k}{3}=[0 5*10^6];
        patient{1}{j,k}{4}=0;
        patient{1}{j,k}{5}=0;
        patient{1}{j,k}{6}='n'; 
        patient{1}{j,k}{7}='n';
        patient{1}{j,k}{8}='rr';
        patient{1}{j,k}{9}=0;
        patient{1}{j,k}{10}=0; 
        patient{1}{j,k}{11}='f';

    end
    m=n6+1;
    for k=m:n7
        %patient{1}{j,k}=newpatientR;
        patient{1}{j,k}=cell(1,11);
        patient{1}{j,k}{1}=0;
        patient{1}{j,k}{2}='r';
        patient{1}{j,k}{3}=[10^7 5*10^6];
        patient{1}{j,k}{4}=0;
        patient{1}{j,k}{5}=0;
        patient{1}{j,k}{6}='n'; 
        patient{1}{j,k}{7}='n';
        patient{1}{j,k}{8}='nr';
        patient{1}{j,k}{9}=0;
        patient{1}{j,k}{10}=0; 
        patient{1}{j,k}{11}='f';
    end
    m=n7+1;
    for k=m:nbp
        %patient{1}{j,k}=newpatientU;
        patient{1}{j,k}=cell(1,11);
        patient{1}{j,k}{1}=0;
        patient{1}{j,k}{2}='u';
        patient{1}{j,k}{3}=[0 0];
        patient{1}{j,k}{4}=0;
        patient{1}{j,k}{5}=0;
        patient{1}{j,k}{6}='n'; 
        patient{1}{j,k}{7}='n';
        patient{1}{j,k}{8}='ui';
        patient{1}{j,k}{9}=0;
        patient{1}{j,k}{10}=0; 
        patient{1}{j,k}{11}='f';

    end
end

% We compute the age of infection 

for j=1:nbs
    for k=1:nbp
        if (patient{1}{j,k}{2}=='r' | patient{1}{j,k}{2}=='n')
            age=valueexprnd(expinit);
            expinit=nextvalueexprnd(expinit);
            while (age>0)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % We compute the new dose of the patient after a time step
                % =timestepinfectionage
                %patient{1}{j,k}=evolutiondosepatient(shift,bt,et,TN,TR,tau,kappa,gamma,betaN1,betaR1,betaN2,betaR2,betaN3,betaR3,patient{1}{j,k});
                if (patient{1}{j,k}{2}~='u') 
  
                    a1=patient{1}{j,k}{1};
                    a2=a1+timestep/2;
                    a3=a1+timestepinfectionage;
                    T=[a1 a2 a3];
                    W=patient{1}{j,k}{3};
    
                    V=zeros(1,2);
                    epsilon=timestepinfectionage/100;
                    c1=bt-epsilon;
                    c2=et-epsilon;
    
                    if (a1<c1 | a1>=c2)    % no antibiotic treatment period
                    [t,x]=ode45(@Bacterialload,T,W,[],tau,kappa,gamma,TREC,betaN1,betaR1,betaN3,betaR3);
                    V(1)=x(3,1);
                    V(2)=x(3,2);
                    else % antibiotic treatment period
                    [t,x]=ode45(@Bacterialload,T,W,[],tau,kappa,gamma,TREC,betaN2,betaR2,betaN3,betaR3);
                    V(1)=x(3,1);
                    V(2)=x(3,2);
                    end 
     
                    % we compute the bacterial load after 1 timestep
                    patient{1}{j,k}{3}=V;

                    % we compute the age of infection
                    patient{1}{j,k}{1}=patient{1}{j,k}{1}+timestepinfectionage;
        
                    % We compute the new infectiousness status of the patient
                    if (patient{1}{j,k}{3}(1)>=TN)
                        patient{1}{j,k}{6}='i';
                    else 
                        patient{1}{j,k}{6}='n';  
                    end
                    if (patient{1}{j,k}{3}(2)>=TR)
                        patient{1}{j,k}{7}='i';
                    else 
                        patient{1}{j,k}{7}='n';   
                    end   
    
                end 

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %patient{1}{j,k}=evolutionpatient(timestepinfectionage,patient{1}{j,k});
                patient{1}{j,k}{4}=patient{1}{j,k}{4}-timestepinfectionage;
                patient{1}{j,k}{5}=patient{1}{j,k}{5}+timestepinfectionage;
                age=age-timestepinfectionage;
            end
        end
    end
end



% we compute the remaining the stay in the hospital  

for j=1:nbs

        for k=1:nbp

            if (patient{1}{j,k}{2}=='u')
                patient{1}{j,k}{4}=valueexprnd(expexitU);       
                expexitU=nextvalueexprnd(expexitU); 
            end 
            if (patient{1}{j,k}{2}=='n')
                patient{1}{j,k}{4}=valueexprnd(expexitN);
                expexitN=nextvalueexprnd(expexitN);
            end 
            if (patient{1}{j,k}{2}=='r')
                patient{1}{j,k}{4}=valueexprnd(expexitR);
                expexitR=nextvalueexprnd(expexitR);
            end 


        end

end 


% patientindexinit will be use the determine if the list of patient free of hcw 
for k=1:nbp
    indexpatientfreeinit(k)=k;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     MAIN PART OF THE PROGRAM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nbday
    ptt=patient{i};
    for s=1:nbshift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Begining of shift s        
%       We create new hcws and we randomly reorganize the index of hcw 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %hcw=newpopulationhcw(nbs,nbh);
        for j=1:nbs
            for k=1:nbh
                %hcw{j,k}=newhcw;
                hcw{j,k}=cell(1,6);
                hcw{j,k}{1}='u';
                hcw{j,k}{2}='u';
                hcw{j,k}{3}=0;
                hcw{j,k}{4}=0;
                hcw{j,k}{5}=0;
                hcw{j,k}{6}='u';
            end 
        end
        for t=1:nbtimestep
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % We determine the index of patients which are not visited by some hcw
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j=1:nbs  
                indexpatientfree=indexpatientfreeinit; 
                % we determine the list of patient free
                for k=1:nbh 
                    if (hcw{j,k}{5}~=0)
                        n1=hcw{j,k}{5};
                        n2=1;
                        while (indexpatientfree(n2)~=n1) 
                            n2=n2+1;
                        end
                        indexpatientfree(n2)=[];      
                    end 
                end 
                
                for k=1:nbh
                    % visits         
                    if (hcw{j,k}{5}==0) 
                        nbpv=size(indexpatientfree);
                        aux=randperm(nbpv(2));           
                        n=indexpatientfree(aux(1));
                        indexpatientfree(aux(1))=[];
                        ptt{j,n}{10}=k; 
                        ptt{j,n}{11}=hcw{j,k}{6};   
                        % we compute the duration of the visit
                        visitlength=valueexprnd(expvisit);          
                        expvisit=nextvalueexprnd(expvisit);
                        hcw{j,k}{3}=visitlength;
                        hcw{j,k}{5}=n;  
                    end 
                    % contamination and infection 
                    if (hcw{j,k}{3}<=timestep) % If the visit of the hcw is finished 
                        n=hcw{j,k}{5};  % we compute the index of patient 
                        if (ptt{j,n}{6}=='i'  | ptt{j,n}{7}=='i')
                            x1=valuerandom(random1);   
                            random1=nextvaluerandom(random1);
                                              
                            % contamination part 
                            if (x1<PC) % if there is contamination 
                             % if the patient is infectious  
                             
                                %hcwaux=contamination(hcw{j,k},ptt{j,n});
                                hcwaux=hcw{j,k};
                                if (ptt{j,n}{6}=='i')
                                    hcwaux{1}='c';
                                    hcwaux{6}='n';
                                end
                                if (ptt{j,n}{7}=='i')
                                    hcwaux{2}='c';
                                    hcwaux{6}='r';
                                end

                                % we compute the time of contamination 
                                tc=valueexprnd(expcontamhcw);    
                                expcontamhcw=nextvalueexprnd(expcontamhcw);
                                hcwaux{4}=tc;
                        
                            else
                                hcwaux=hcw{j,k};
                            end
                        else 
                            hcwaux=hcw{j,k};
                        end
                
                    
                        % infection part 
                        if (hcw{j,k}{6}~='u' & ptt{j,n}{2}~='r')
                            x2=valuerandom(random2);   
                            random2=nextvaluerandom(random2);                                  
                            % infection part 
                            if (x2<PI) % if there is infection                                                                               
                                charaux=ptt{j,n}{2};

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %ptt{j,n}=infection(hcw{j,k},ptt{j,n});
                                if (hcw{j,k}{1}=='c' & hcw{j,k}{2}=='u')
                                    if (ptt{j,n}{2}=='u')
                                        ptt{j,n}{2}='n';
                                        ptt{j,n}{3}=[10^7 0];
                                    end
                                end
                                if (hcw{j,k}{1}=='u' & hcw{j,k}{2}=='c')
                                    if (ptt{j,n}{2}=='u' )
                                        ptt{j,n}{2}='r';
                                        ptt{j,n}{3}=[0 5*10^6];
                                        ptt{j,n}{8}='rr';
                                    end
                                    if (ptt{j,n}{2}=='n')
                                        ptt{j,n}{2}='r';
                                        ptt{j,n}{3}=[10^10 5*10^6];
                                        ptt{j,n}{8}='rs';
                                    end
                                end
                                if (hcw{j,k}{1}=='c' & hcw{j,k}{2}=='c')
                                    if (ptt{j,n}{2}=='u' )
                                        ptt{j,n}{2}='r';
                                        ptt{j,n}{3}=[10^7 5*10^6];
                                        ptt{j,n}{8}='nr';
                                    end
                                    if (ptt{j,n}{2}=='n')
                                        ptt{j,n}{2}='r';
                                        ptt{j,n}{3}=[10^10 5*10^6];
                                        ptt{j,n}{8}='rs';
                                    end
                                end
   
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                
                                
                            
                                % we update the time remaining to stay in the hospital    
                                if (charaux=='u' & ptt{j,n}{2}=='n')
                                    ptt{j,n}{4}=valueexprnd(expexitN);
                                    expexitN=nextvalueexprnd(expexitN);
                                end
                                if (charaux=='u' & ptt{j,n}{2}=='r')
                                    ptt{j,n}{4}=valueexprnd(expexitR);
                                    expexitR=nextvalueexprnd(expexitR);
                                end  
                                if (charaux=='n' & ptt{j,n}{2}=='r')
                                    ptt{j,n}{4}=valueexprnd(expexitR);
                                    expexitR=nextvalueexprnd(expexitR);
                                end  
                              
                            
                            end %if (x2<PI)
                        end %if (hcw{6}~='u' & patient{j,n}{2}~='r')
                        
                        
                        hcw{j,k}=hcwaux;
                        ptt{j,n}{10}=0; % the patient stop to visit the hcw 
                        ptt{j,n}{11}='f';
                            
                                   
                    end % if (hcw{j,k}{3}<=timestep)
                end %for k=1:nbh
            end %for j=1:nbs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Evolution of patients and hcws at the end of time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            for j=1:nbs
                for k=1:nbp
                    %ptt{j,k}=evolutionpatient(timestep,ptt{j,k});
                    ptt{j,k}{4}=ptt{j,k}{4}-timestep;
                    ptt{j,k}{5}=ptt{j,k}{5}+timestep;
                    
                    if (ptt{j,k}{4}<=0)     % if the patient leave the hospital 
                        n=ptt{j,k}{10};
                        if (n>0) % if the patient is visited by some hcw 
                            hcw{j,n}{3}=0; 
                            hcw{j,n}{5}=0;
                        end
                        %ptt{j,k}=newpatientU;
                        ptt{j,k}=cell(1,11);
                        ptt{j,k}{1}=0;
                        ptt{j,k}{2}='u';
                        ptt{j,k}{3}=[0 0];
                        ptt{j,k}{4}=0;
                        ptt{j,k}{5}=0;
                        ptt{j,k}{6}='n'; 
                        ptt{j,k}{7}='n';
                        ptt{j,k}{8}='ui';
                        ptt{j,k}{9}=0;
                        ptt{j,k}{10}=0; 
                        ptt{j,k}{11}='f';

                        ptt{j,k}{4}=valueexprnd(expexitU);
                        expexitU=nextvalueexprnd(expexitU);
                        
            
                    end
                end
            end  
            for j=1:nbs
                for k=1:nbh
                    %hcw{j,k}=evolutionhcw(timestep,hcw{j,k});
                    hcw{j,k}{3}=hcw{j,k}{3}-timestep;
                    hcw{j,k}{4}=hcw{j,k}{4}-timestep;


                    if (hcw{j,k}{3}<=0) % if the time remaining for the visit is negative then it is egal to 0
                        hcw{j,k}{3}=0;  
                        hcw{j,k}{5}=0;
                    end
                    if (hcw{j,k}{4}<=0)  % if the time remaining to be contaminated is <=0  
                                %then the hcw becomes uncontaminated
                        hcw{j,k}{4}=0;    
                        hcw{j,k}{1}='u';
                        hcw{j,k}{2}='u';
                        hcw{j,k}{6}='u';

                    end


                end
            end
            
            
            
            
            
            
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       End of the time steps 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        end %for t=1:nbtimestep
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Evolution of the infectiousness of patients at the end of the shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

        for j=1:nbs
            for k=1:nbp
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %ptt{j,k}=evolutiondosepatient(shift,bt,et,TN,TR,tau,kappa,gamma,betaN1,betaR1,betaN2,betaR2,betaN3,betaR3,ptt{j,k});
                if (ptt{j,k}{2}~='u') 
  
                    a1=ptt{j,k}{1};
                    a2=a1+timestep/2;
                    a3=a1+timestepinfectionage;
                    T=[a1 a2 a3];
                    W=ptt{j,k}{3};
    
                    V=zeros(1,2);
                    epsilon=timestepinfectionage/100;
                    c1=bt-epsilon;
                    c2=et-epsilon;
    
                    if (a1<c1 | a1>=c2)    % no antibiotic treatment period
                    [t,x]=ode45(@Bacterialload,T,W,[],tau,kappa,gamma,TREC,betaN1,betaR1,betaN3,betaR3);
                    V(1)=x(3,1);
                    V(2)=x(3,2);
                    else % antibiotic treatment period
                    [t,x]=ode45(@Bacterialload,T,W,[],tau,kappa,gamma,TREC,betaN2,betaR2,betaN3,betaR3);
                    V(1)=x(3,1);
                    V(2)=x(3,2);
                    end 
     
                    % we compute the bacterial load after 1 timestep
                    ptt{j,k}{3}=V;

                    % we compute the age of infection
                    ptt{j,k}{1}=ptt{j,k}{1}+timestepinfectionage;
        
                    % We compute the new infectiousness status of the patient
                    if (ptt{j,k}{3}(1)>=TN)
                        ptt{j,k}{6}='i';
                    else 
                        ptt{j,k}{6}='n';  
                    end
                    if (ptt{j,k}{3}(2)>=TR)
                        ptt{j,k}{7}='i';
                    else 
                        ptt{j,k}{7}='n';   
                    end   
    
                end 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end
        end  
            

    
         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       End of Shift s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    end 

 
    patient{i+1}=ptt;
       

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       End of the day i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               WE NOW PLOT THE TRAJECTORIES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nbday=nbday+1;
result=zeros(nbday,8);    % historypatient contains the number of patient per 
                                            %day in each dialyze center and for each statue S, IN, 
                                            % IR, and IS. 
time=zeros(nbday,1);
t=zeros(nbday,1);
for i=1:nbday
    time(i,1)=i-1;
    for j=1:nbs
        for k=1:nbp
            if (patient{i}{j,k}{2}=='u')
                result(i,1)=result(i,1)+1;
            end
            if (patient{i}{j,k}{2}=='n')
                result(i,2)=result(i,2)+1;
            end
            if (patient{i}{j,k}{2}=='r')
                result(i,3)=result(i,3)+1;            end
            if (patient{i}{j,k}{8}=='nr')
               result(i,4)=result(i,4)+1;
            end   
            if (patient{i}{j,k}{8}=='rs')
               result(i,5)=result(i,5)+1;
            end  
            if (patient{i}{j,k}{8}=='rr')
               result(i,6)=result(i,6)+1;
            end 
            if (patient{i}{j,k}{6}=='i')
               result(i,7)=result(i,7)+1;
            end
            if (patient{i}{j,k}{7}=='i')
               result(i,8)=result(i,8)+1;
            end
                
        end
    end
end

nbptotal=nbs*nbp;                                                       
result=100*result/nbptotal; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               WE PLOT THE S(t), IN(t), and IR(t)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result(:,1)=result(:,2)+result(:,3);
h=plot(time,result(:,1),time,result(:,2),time,result(:,3),time,result(:,4),time,result(:,5),time,result(:,6),'LineWidth',1.5);axis([-Inf Inf 0 100])
h=legend('% of N and R','% of N','% of R','% of NR','% of RS','% of RR');
xlabel('days');
ylabel('Percentage of patients');
AV=AV*(60*24);
title(['Nbs=',num2str(nbs),', AV=',num2str(AV),' mn, bt=',num2str(bt),' days, et=',num2str(et),' days']);
axis([-Inf Inf 0 100])
saveas(h,filename)
NCLOCK=fix(clock)