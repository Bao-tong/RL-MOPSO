function Position = Leaderparticle(pop,i,nPop)

f1=pop(i).Explore(1);
f2=i;
f3=pop(i).Explore(2);

    IPosition=LeaderFindPosition1(pop,f1,f2,f3);
    Position=IPosition;




end