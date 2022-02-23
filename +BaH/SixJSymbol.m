function[res]=SixJSymbol(j1,j2,j3,j4,j5,j6)

    import TlF.ThreeJSymbol

    if mod(2*j1,1)~=0 || mod(2*j2,1)~=0 || mod(2*j3,1)~=0 || mod(2*j4,1)~=0 || mod(2*j5,1)~=0 || mod(2*j6,1)~=0
		res=0;
        return
    end

	m1=-j1:1:j1;
    m2=-j2:1:j2;
    m3=-j3:1:j3;
    m4=-j4:1:j4;
    m5=-j5:1:j5;
    m6=-j6:1:j6;


	res=0;

	for M1=m1
		for M2=m2
			for M3=m3
				for M4=m4
					for M5=m5
						for M6=m6
							power=j1+j2+j3+j4+j5+j6-M1-M2-M3-M4-M5-M6;
							res=res+((-1)^power)*ThreeJSymbol(j1,-M1,j2,-M2,j3,-M3)*ThreeJSymbol(j1,M1,j5,-M5,j6,M6)*ThreeJSymbol(j4,M4,j2,M2,j6,-M6)*ThreeJSymbol(j4,-M4,j5,M5,j3,M3);
                        end
                    end
                end
            end
        end
    end
	
end