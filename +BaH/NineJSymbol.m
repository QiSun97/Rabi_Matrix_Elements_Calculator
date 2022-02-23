function [res]=NineJSymbol(j1,j2,j3,j4,j5,j6,j7,j8,j9)
    
    import TlF.ThreeJSymbol

    if mod(2*j1,1)~=0 || mod(2*j2,1)~=0 || mod(2*j3,1)~=0 || mod(2*j4,1)~=0 || mod(2*j5,1)~=0 || mod(2*j6,1)~=0 || mod(2*j7,1)~=0 || mod(2*j8,1)~=0 || mod(2*j9,1)~=0
		res=0;
        return
    end

	m1=-j1:1:j1;
    m2=-j2:1:j2;
    m3=-j3:1:j3;
    m4=-j4:1:j4;
    m5=-j5:1:j5;
    m6=-j6:1:j6;
    m7=-j7:1:j7;
    m8=-j8:1:j8;
    m9=-j9:1:j9;


	res=0;

	for M1=m1
		for M2=m2
			for M3=m3
				for M4=m4
					for M5=m5
						for M6=m6
                            for M7=m7
                                for M8=m8
                                    for M9=m9
                                        res=res+ThreeJSymbol(j1,M1,j2,M2,j3,M3)*ThreeJSymbol(j4,M4,j5,M5,j6,M6)*ThreeJSymbol(j7,M7,j8,M8,j9,M9)*ThreeJSymbol(j1,M1,j4,M4,j7,M7)*ThreeJSymbol(j2,M2,j5,M5,j8,M8)*ThreeJSymbol(j3,M3,j6,M6,j9,M9);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
   
end