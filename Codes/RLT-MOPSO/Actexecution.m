function pop=Actexecution(rt,pop,i,VarSize,c1,c2,VarMin,VarMax,it,nPop,rows,cols,w,rep,beta,nVar)
CostFunction = @(x) MOP2(x);            
if rt==1
       pop(i).Leader.Position = Leaderparticle(pop,i,nPop,rows, cols);
       pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(pop(i).Leader.Position-pop(i).Position);
        pop(i).topology=1;
        end
        
        if rt==2
       pop(i).Best.Position = Memoryparticle(pop,i,nPop);
       pop(i).Leader.Position = Leaderparticle1(pop,i,nPop);       
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(pop(i).Leader.Position-pop(i).Position);
        pop(i).topology=2;
        end
        
         if rt==3
        leader = SelectLeader(rep, beta);       
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);      
        pop(i).topology=3;
         end        
     
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        


for p=1:nVar
        pop(i).Position(p) = max(pop(i).Position(p), VarMin(p));
        pop(i).Position(p) = min(pop(i).Position(p), VarMax(p));
end
        
        pop(i).Cost = CostFunction(pop(i).Position);
end