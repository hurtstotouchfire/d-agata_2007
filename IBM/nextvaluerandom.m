function z=nextvaluerandom(y)

z=cell(1,2);
n=size(y{2});
if (y{1}<n(1))
    y{1}=y{1}+1;
    z=y;
else
    z=newrandom;
end

    