function I=evolutionpatient(timestep,bt,et,TN,TR,TREC,tau,kappa,gamma,betaN1,betaR1,betaN2,betaR2,betaN3,betaR3,patient)



if (patient{2}~='u') 
  
    a1=patient{1};
    a2=a1+timestep/2;
    a3=a1+timestep;
    T=[a1 a2 a3];
    W=patient{3};
    
    V=zeros(1,2);
    epsilon=timestep/100;
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
    patient{3}=V;

    % we compute the age of infection
    patient{1}=patient{1}+timestep;
        
    % We compute the new infectiousness status of the patient
    if (patient{3}(1)>=TN)
        patient{6}='i';
    else 
        patient{6}='n';  
    end
    if (patient{3}(2)>=TR)
        patient{7}='i';
    else 
        patient{7}='n';   
    end   
    
end 


I=patient;