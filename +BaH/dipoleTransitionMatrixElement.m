function[x]=dipoleTransitionMatrixElement(G,E,p)
    
    import BaH.*
    
    if G(1)==0 || G(1)==2
        [G,Cg]=project(G);
    elseif G(1)==1
        L=G(2);
        S=G(3);
        sig=G(4);
        W=G(5);
        J=G(6);
        F=G(7);
        MF=G(8);
        P=G(9);
        if L-sig==W
            G=[[L,S,-sig,W,J,F,MF,P];[-L,S,sig,-W,J,F,MF,P]];
        else
            G=[[L,S,sig,W,J,F,MF,P];[-L,S,-sig,-W,J,F,MF,P]];
        end

        if P==1
            Cg=[1/sqrt(2),1/sqrt(2)*(-1)^(J-S)];
        else
            Cg=[1/sqrt(2),-1/sqrt(2)*(-1)^(J-S)];
        end
    end
    
    if E(1)==0 || E(1)==2
        [E,Ce]=project(E);
    elseif E(1)==1
        L=E(2);
        S=E(3);
        sig=E(4);
        W=E(5);
        J=E(6);
        F=E(7);
        MF=E(8);
        P=E(9);
        if L-sig==W
            E=[[L,S,-sig,W,J,F,MF,P];[-L,S,sig,-W,J,F,MF,P]];
        else
            E=[[L,S,sig,W,J,F,MF,P];[-L,S,-sig,-W,J,F,MF,P]];
        end

        if P==1
            Ce=[1/sqrt(2),1/sqrt(2)*(-1)^(J-S)];
        else
            Ce=[1/sqrt(2),-1/sqrt(2)*(-1)^(J-S)];
        end
    end
    
    
    
    x=0;
    for i=1:2
        for j=1:2
            y=0;
            g=G(i,:);
            e=E(j,:);
            
            
            Le=e(1);
            Se=e(2);
            sige=e(3);
            We=e(4);
            Je=e(5);
            Fe=e(6);
            MFe=e(7);
            Pe=e(8);
 
            Lg=g(1);
            Sg=g(2);
            sigg=g(3);
            Wg=g(4);
            Jg=g(5);
            Fg=g(6);
            MFg=g(7);
            Pg=g(8);
            
            
            power=Fe-MFe+Fg+Jg+1/2+1;
            y=(-1)^power*sqrt((2*Fe+1)*(2*Fg+1))*ThreeJSymbol(Fe,-MFe,1,p,Fg,MFg)*SixJSymbol(Jg,Fg,1/2,Fe,Je,1)*(1-delta_kr(Pe,Pg));
            z=0;
            for q=[-1,0,1]
                z=z+(-1)^(Je-We)*sqrt((2*Je+1)*(2*Jg+1))*ThreeJSymbol(Je,-We,1,q,Jg,Wg)*delta_kr(sigg,sige);
            end

            y=y*z;
            x=x+y*Cg(i)*Ce(j);
        end
    end

end