function dy=EDOpatientinfection1(t,y,NUN,NUR,NUV,BETAV,PI,HNN,HNR,HRR)

%
% bn=nbh/nbp*Hn
% bu=nbh/nbp*Hu
% bf=1-(nbh/nbp)
% 
%

dy=zeros(5,1);  
dy(1)=NUN*y(2)+NUR*(y(3)+y(4)+y(5))-NUV*BETAV*PI*(HNN+HNR+HRR)*y(1);
dy(2)=NUV*BETAV*PI*HNN*y(1)-NUV*BETAV*PI*(HNR+HRR)*y(2)-NUN*y(2);
dy(3)=NUV*BETAV*PI*(HNR+HRR)*y(2)-NUR*y(3);
dy(4)=NUV*BETAV*PI*HRR*y(1)-NUR*y(4);
dy(5)=NUV*BETAV*PI*HNR*y(1)-NUR*y(5);
