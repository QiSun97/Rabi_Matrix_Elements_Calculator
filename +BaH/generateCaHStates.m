function[XE,BE,nground,nexcite]=generateCaHStates()
%Ns: XSigma states rotational quantum number N
%Omega: angular momentum projection for the APi state with N=0
%Parity:
%NsB: BSigma states rotational quantum number N
%Outputs: [X/A/B,N,S,J,F,Mf] for Hund case (b)
%[X/A/B,lambda,S,sigma,omega,J,F,Mf,parity] for Hund case (A)
%BaH is assumed to have bosonic Ba isotope with I=0

XE=[];
BE=[];
nground = 0;
nexcite = 0;
%XSigma states: S=1/2, L=0 (Hund's case b)
% here i only use N=1, J=3/2, F=2
n=1;
s=1/2;
for j = [1/2, 3/2]
    for f=[j-1/2, j+1/2]
        for m=-f:1:f
            XE=[XE;[0,n,s,j,f,m]]; %[X/A/B,N,S,J,F,Mf]
            nground = nground + 1;
        end
    end
end

for j=[1/2]
    for i=[-1/2,1/2]
        f=j+i;
        for m=-f:1:f
            BE=[BE;[1,1,1/2,1/2,1/2,j,f,m,1]];%[X/A/B,lambda,S,sigma,omega,J,F,Mf,parity] (A^2Pi_1/2 state -> L=1, S=1/2, )
            nexcite = nexcite + 1;
        end
    end
end

%BSigma states: S=1/2, L=0 (Hund's case b)
% here i only use N=0, J=1/2, F=1
% n=0;
% s=1/2;
% j=1/2;
% for f=[j+1/2]
%     for m=-f:1:f
%         BE=[BE;[2,n,s,j,f,m]]; %[X/A/B,N,S,J,F,Mf]
%         nexcite = nexcite + 1;
%     end
% end
end