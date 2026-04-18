function C=IR(rep,rop,Obj)

a=numel(rep);
for i=1:a
    for j=1:Obj
PopObj_current(i,j)=rep(i).Cost(j);
    end
end

a=numel(rop);
for i=1:a
    for j=1:Obj
PopObj_previous(i,j)=rop(i).Cost(j);
    end
end

Distance = pdist2(PopObj_current, PopObj_previous);

minDistCurrent = min(Distance, [], 2);

minDistPrevious = min(Distance, [], 1);

ImprovementRate = mean(minDistPrevious) - mean(minDistCurrent);
C=ImprovementRate;


