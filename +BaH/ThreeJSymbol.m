function[res]=ThreeJSymbol(j1,m1,j2,m2,j3,m3)
    import TlF.cg
    res=((-1)^(j1-j2-m3)/sqrt(2*j3+1))*cg(j1,m1,j2,m2,j3,-m3);
end