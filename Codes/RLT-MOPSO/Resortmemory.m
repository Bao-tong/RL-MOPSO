function pop=Resortmemory(pop,nPop,Obj, EPnum)


for e=1:EPnum

for j=1:Obj
for i=1:nPop
    A(i,j)=pop(i).Best.Cost(j);
end
end

B(1,:)=max(A);
B(2,:)=min(A);

for i=1:nPop
    for j=1:Obj
        N(j)=(pop(i).Best.Cost(j)-B(2,j))/(B(1,j)-B(2,j)); 
        if j==1
            N=[N(1)];
            if j>1
            N=[N,N(j)];   
            end
        end
    end
             pop(i).Normalization=N;
end

if e==1
for i=1:nPop
     for u=1:nPop
        for j=1:Obj
    C(i,u)=sum(sqrt(pop(i).Best.Cost(j)-pop(u).Best.Cost(j)));      
        end
     end  
end
end

min_nonzero_vals = nan(size(C, 1), 1);


for row = 1:size(C, 1)
    non_zero_elements = C(row, C(row, :) ~= 0);
    
    if ~isempty(non_zero_elements)
        min_nonzero_vals(row) = min(non_zero_elements);
    end
end
for i=1:nPop
        D(i,e)=find(min_nonzero_vals(i)==C(i,:),1,'first');
end

for i=1:nPop
    C(i,D(i))=inf;
end 

if e==1
    for i=1:nPop
    pop(i).Memory=D(i,1);
    end
end

if e==2
    for i=1:nPop
        pop(i).Memory=[D(i,1),D(i,2)];
    end
end
end