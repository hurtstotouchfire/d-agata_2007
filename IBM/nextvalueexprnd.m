function z=nextvalueexprnd(y)

z=cell(1,3);
n=size(y{2});
if (y{1}<n(1))
    y{1}=y{1}+1;
    z=y;
else
    z=newexprnd(y{3});
    
end
