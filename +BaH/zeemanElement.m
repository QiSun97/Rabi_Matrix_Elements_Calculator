function[x]=zeemanElement(G,E,p) %they both have to be in case (b); L is assumed to be 0 (Sigma state)

    import BaH.*
    if G(1)~=E(1)
        x=0;
        return
    end
    if G(1)==1
        x=0;
        return
    else
        Ng=G(2);
        Sg=G(3);
        Jg=G(4);
        Fg=G(5);
        MFg=G(6);
        Ne=E(2);
        Se=E(3);
        Je=E(4);
        Fe=E(5);
        MFe=E(6);
        power=Fe-MFe+Fg+Je+1+1/2+Je+Ne+1+Se;

        x=(-1)^power*sqrt((2*Fe+1)*(2*Fg+1))*sqrt((2*Je+1)*(2*Jg+1))*sqrt(Se*(Se+1)*(2*Se+1))*delta_kr(Se,Sg)*delta_kr(Ne,Ng)*ThreeJSymbol(Fe,-MFe,1,p,Fg,MFg)*SixJSymbol(Fe,Je,1/2,Jg,Fg,1)*SixJSymbol(Je,Se,Ne,Se,Jg,1);
    end
end