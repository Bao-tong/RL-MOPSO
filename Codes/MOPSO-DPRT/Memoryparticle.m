function Position = Memoryparticle(pop,i,nPop)


f1=pop(i).Memory(1);
f2=i;
f3=pop(i).Memory(2);

IPosition=MemoryFindPosition(pop,f1,f2,f3);
Position=IPosition;




end