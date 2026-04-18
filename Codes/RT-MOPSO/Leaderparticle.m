function Position = Leaderparticle(pop,i,nPop,rows, cols)

pop(i).neighbors=TwoDim(rows, cols, nPop,i);

if numel(pop(i).neighbors)==2

    f1=pop(i).neighbors(1);
    f2=i;
    f3=pop(i).neighbors(2);
   [IPosition,bestparticle]=LeaderFindPosition(pop,f1,f2,f3);
    Position=IPosition;
    
end

if numel(pop(i).neighbors)==3

    f1=pop(i).neighbors(1);
    f2=i;
    f3=pop(i).neighbors(2);
   [IPosition,bestparticle]=LeaderFindPosition(pop,f1,f2,f3);

    A=bestparticle;
    B=pop(i).neighbors(3);
if Dominates(pop(A), pop(B))
    C=1 ;
else
    C=0;
    IPosition=pop(B).Position;
end
    Position=IPosition;
end



if numel(pop(i).neighbors)==4

    f1=pop(i).neighbors(1);
    f2=i;
    f3=pop(i).neighbors(2);
   [IPosition,bestparticle]=LeaderFindPosition(pop,f1,f2,f3);

    f1=pop(i).neighbors(3);
    f2=bestparticle;
    f3=pop(i).neighbors(4);

    
       [IPosition,bestparticle]=LeaderFindPosition(pop,f1,f2,f3);
    Position=IPosition;
    
end

end