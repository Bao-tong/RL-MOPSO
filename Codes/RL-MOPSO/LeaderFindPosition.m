function [IPosition,bestparticle]=LeaderFindPosition(pop,f1,f2,f3)

            if Dominates(pop(f2), pop(f1))
               A=1;
            else
                A=0;
            end
            
            if Dominates(pop(f1), pop(f2))
               B=1;
            else
                B=0;
            end
                 
            if Dominates(pop(f2), pop(f3))
               C=1;
            else
                C=0;
            end
            
            if Dominates(pop(f3), pop(f2))
               D=1;
            else
                D=0;
            end
            
            if Dominates(pop(f1), pop(f3))
               E=1;
            else
                E=0;
            end
            
            if Dominates(pop(f3), pop(f1))
               F=1;
            else
                F=0;
            end
            
           if (A==1&&C==1)||(A==1&&E==1)|| (C==1&&F==1)        
                IPosition=pop(f2).Position;
                bestparticle=f2;
            end
            if (B==1&& E==1)||(B==1&&C==1)|| (D==1&&E==1) 
               IPosition=pop(f1).Position;
               bestparticle=f1;
            end
            if (D==1 && F==1)||(A==1&&D==1)|| (B==1&&F==1)
               IPosition=pop(f3).Position;
               bestparticle=f3;
            end
            
            if (A==1 && F==1&&C==0&&D==0)||(A==1&&C==0&&D==0&&E==0&&F==0)||(F==1&&A==0&&B==0&&C==0&&D==0)
            bestparticle=randsample([f2,f3], 1);
                IPosition=pop(bestparticle).Position;
            end
            if (C==1 && E==1&&A==0&&B==0)||(C==1&&A==0&&B==0&&E==0&&F==0)||(E==1&&A==0&&B==0&&C==0&&D==0)
             bestparticle=randsample([f2,f1], 1);
                IPosition=pop(bestparticle).Position;
            end
            if (B==1 && D==1&&E==0&&F==0)||(B==1&&C==0&&D==0&&E==0&&F==0)||(D==1&&A==0&&B==0&&E==0&&F==0)
              bestparticle=randsample([f3,f1], 1);
                IPosition=pop(bestparticle).Position;
            end            
             if A==0 && B==0&&C==0&&D==0&&E==0&&F==0
             bestparticle=randsample([f2,f3,f1], 1);
                 IPosition=pop(bestparticle).Position;
            end           


end