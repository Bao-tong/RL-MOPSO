function pop=repdis(pop,nPop,rep,Obj)


for j=1:Obj
for i=1:nPop
    A(i,j)=pop(i).Cost(j);
end
end

B(1,:)=max(A);
B(2,:)=min(A);

for i=1:nPop
    for j=1:Obj
        N(j)=(pop(i).Cost(j)-B(2,j))/(B(1,j)-B(2,j)); 
        if j==1
            N=[N(1)];
            if j>1
            N=[N,N(j)];   
            end
        end
    end
      pop(i).Normalization=N;
end


for i=1:nPop
     for u=1:numel(rep)
        for j=1:Obj
    C(i,u)=sum(sqrt((pop(i).Cost(j))^2-(rep(u).Cost(j))^2));      
        end
     end  
end


numRows=nPop;
numCols=numel(rep);
result = zeros(numRows, 1);

for i = 1:numRows
    if all(C(i, :) == 0)
        result(i) = 1;
    else
        nonZeroValues = C(i, C(i, :) ~= 0);
        [minValue, minIdx] = min(nonZeroValues);
        result(i) = find(C(i, :) == minValue, 1, 'first');
    end
    pop(i).repdis= result(i);
end
