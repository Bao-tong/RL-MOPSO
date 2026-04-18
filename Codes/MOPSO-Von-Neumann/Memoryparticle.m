function Position = Memoryparticle(pop,i,nPop,rows, cols)


pop(i).neighbors=TwoDim(rows, cols, nPop,i);






if i==1
    f1=nPop;
    f2=i;
    f3=i+1;

    IPosition=MemoryFindPosition(pop,f1,f2,f3);
    Position=IPosition;
end

if i>1&&i<nPop
    f1=i-1;
    f2=i;
    f3=i+1;
    IPosition=MemoryFindPosition(pop,f1,f2,f3);
    Position=IPosition;
end

if i==nPop
    f1=i-1;
    f2=i;
    f3=1;
    IPosition=MemoryFindPosition(pop,f1,f2,f3);
    Position=IPosition;
end


end