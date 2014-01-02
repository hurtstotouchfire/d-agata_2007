function I=newpatientR
%
% This function creates a new suceptible individual. 
%
%
% I{1}= a the infection age of the inidividual  
%
% I{2}= 'u' (uninfected), 'n' (infected only by non resistant strain),
%        or 'r' (infected by both resistant and non resistant strain)
%
% I{3}=[V-(a) V+(a)] the bacterial load
% I{4}= remaining to time spend in the hospital 
% I{5}= time already spend in the hospital 
%
% I{6}=    'n' (non infectious for non resistant strain), 'i' (infectious for non resistant strain) 
% 
% I{7}=    'n' (non infectious for resistant strain), 'i' (infectious for resistant strain) 
%
% I{8}=    'ui' uninfected by resistant strain, 
%          'nr' the ones that were uninfected and were infected by both resistant non-resistant 
%           strain,  
%
%          'rs' the ones that were infected by non-resistant strain, and then infected 
%           by resistant strain, 
%
%          'rr' the ones that were uninfected and were infected only by resistant strain.  
%     
%     
% I{9}=    time remaining for a visit of an heathcare workers
%          Remark: this quatity is <=0 if there is no visit
%
%
% I{10}=    is the index of the hcw if the patient is visited, and 0 otherwise
%        
%
% I{11}=   is the type of hcw visited 
%          'f'  if the patient is hcw free 
%          'u' if the patient is visited by an hcw of type 'u'
%          'n' if the patient is visited by an hcw of type 'n'
%          'r' if the patient is visited by an hcw of type 'r'
%     
I=cell(1,11);
I{1}=0;
I{2}='r';
I{3}=[10^10 5*10^6];
I{4}=0;
I{5}=0;
I{6}='n'; 
I{7}='n';
I{8}='rs';
I{9}=0;
I{10}=0; 
I{11}='f';





