function rt=RL(pop,rep,nPop,i,it,MaxIt)

greedy=0.5/(1+exp((10*(it-0.6*MaxIt))/MaxIt));
a=rand;
if a>1-greedy
rt=randi([1, 3]);
else
    A(1,1)=Spacing(rep);
    A(2,1)=IR(rep);   
    rep = [rep
        pop(i)
         pop(~[pop.IsDominated])];
    
    rep = DetermineDomination(rep);
    rep = rep(~[rep.IsDominated]);   
    Grid = CreateGrid(rep, nGrid, alpha);
    for i = 1:numel(rep)
        rep(i) = FindGridIndex(rep(i), Grid);
    end    
    if numel(rep)>nRep
        Extra = numel(rep)-nRep;
        for e = 1:Extra
            rep = DeleteOneRepMemebr(rep, gamma);
        end        
    end
    A(1,2)=Spacing(rep);
    A(2,2)=IR(rep);   
Act=pop(i).topology;
state=pop(i).state;
    if A(1,1)<A(1,2)&&A(2,1)>A(2,2)
        sta=1;
        r=-0.1;
    else if A(1,1)>=A(1,2)&&A(2,1)<=A(2,2)
        sta=2;     
         r=0.1;
        else
        sta=3;
         r=0;
        end
    end
rowValues = pop(i).Q(sta, :);
[~, colIndex] = max(rowValues);
str=num2str(colIndex);
pop(i).Q(state,Act)=(1-0.75)*pop(i).Q(state,Act)+0.75*(r+0.2*pop(i).Q(sta,str)-pop(i).Q(state,Act));
pop(i).state=sta;

rowValues = pop(i).Q(sta, :);
[~, colIndex] = max(rowValues);
str=num2str(colIndex);


if str==1
    rt=1;
end

if str==2
    rt=2;
end

if str==3
    rt=3;
end
end
end