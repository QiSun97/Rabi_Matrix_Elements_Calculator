function[final,coeffs]=project(initial) %Projection of (b) case onto (a) assuming L=0
    import BaH.*

    final=[];
    coeffs=[];
    
    N=initial(2);
    S=initial(3);
    J=initial(4);
    F=initial(5);
    MF=initial(6);
    
    for sig=[-S,S]
        w=sig;
        final=[final;[0,S,sig,w,J,F,MF,(-1)^N]];%[lambda,S,sigma,omega,J,F,Mf] (A^2Pi_1/2 state -> L=1, S=1/2, sigma=+-1/2, lambda=+-1, |omega|=1/2,parity=(-1)^N)
        coeffs=[coeffs,(-1)^(J+w)*sqrt(2*N+1)*ThreeJSymbol(S,sig,N,0,J,-w)];
    end

end