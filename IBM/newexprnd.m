function z=newexprnd(mean)
z=cell(1,3);
z{1}=1;
z{2}=exprnd(mean,1000000,1);
z{3}=mean;