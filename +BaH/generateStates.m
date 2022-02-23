function[XE,AP,BE]=generateStates(Ns,Omega,Parity,NsB)
%Ns: XSigma states rotational quantum number N
%Omega: angular momentum projection for the APi state with N=0
%Parity:
%NsB: BSigma states rotational quantum number N
%Outputs: [X/A/B,N,S,J,F,Mf] for Hund case (b)
%[X/A/B,lambda,S,sigma,omega,J,F,Mf,parity] for Hund case (A)
%BaH is assumed to have bosonic Ba isotope with I=0

if Omega~=1/2 && Omega~=3/2
    disp('Please provide angular momentum projection for the APi state with N=0. It can be either 1/2 or 3/2.')
    return
end


XE=[];
AP=[];
BE=[];

%XSigma states: S=1/2, L=0 (Hund's case b)
for n=Ns
    if n==0
        j=1/2;
        for i=[-1/2,1/2]
            f=j+i;
            for m=-f:1:f
                XE=[XE;[0,n,1/2,j,f,m]]; %[X/A/B,N,S,J,F,Mf]
            end
        end
    else
        for s=[-1/2,1/2]
            j=n+s;
            for i=[-1/2,1/2]
                f=j+i;
                for m=-f:1:f
                    XE=[XE;[0,n,1/2,j,f,m]]; %[X/A/B,N,S,J,F,Mf]
                end
            end
        end
    end
end

%APi states: S=1/2, L=1, sigma=+-1/2, lambda=+-1, |omega|=Omega (1/2 or 3/2), R=0 (Hund's case a),parity=Parity (1 or -1)
for j=Omega
    for i=[-1/2,1/2]
        f=j+i;
        for m=-f:1:f
            AP=[AP;[1,1,1/2,1/2,Omega,j,f,m,Parity]];%[X/A/B,lambda,S,sigma,omega,J,F,Mf,parity] (A^2Pi_1/2 state -> L=1, S=1/2, )
        end
    end
end

%BSigma states: S=1/2, L=0 (Hund's case b)
for n=NsB
    if n==0
        j=1/2;
        for i=[-1/2,1/2]
            f=j+i;
            for m=-f:1:f
                BE=[BE;[2,n,1/2,j,f,m]]; %[X/A/B,N,S,J,F,Mf]
            end
        end
    else
        for s=[-1/2,1/2]
            j=n+s;
            for i=[-1/2,1/2]
                f=j+i;
                for m=-f:1:f
                    BE=[BE;[2,n,1/2,j,f,m]]; %[X/A/B,N,S,J,F,Mf]
                end
            end
        end
    end
end
end