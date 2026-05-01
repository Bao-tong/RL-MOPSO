function D=Spacing(pop,Obj)
a=numel(pop);
for i=1:a
    for j=1:Obj
PopObj(i,j)=pop(i).Cost(j);
    end
end

Distance = pdist2(PopObj, PopObj);

Distance(logical(eye(size(Distance, 1)))) = inf;

d_i = min(Distance, [], 2);

d_mean = mean(d_i);

S = sqrt(sum((d_i - d_mean).^2) / (length(d_i) - 1));
D=S;
end
