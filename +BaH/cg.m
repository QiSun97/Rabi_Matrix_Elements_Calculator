function[res]=cg(j1,m1,j2,m2,J,M)
    assert(isscalar(j1) && isscalar(m1) && isscalar(j2) && isscalar(m2) && isscalar(J) && isscalar(M));
    if j1<0 || j2<0 || J<0 || abs(M)>J || abs(m1)>j1 || abs(m2)>j2 || mod(2*j1,1)~=0 || mod(2*j2,1)~=0 || mod(2*J,1)~=0 || mod(2*m1,1)~=0 || mod(2*m2,1)~=0 || mod(2*M,1)~=0 || j1+m1<0 || j2+m2<0 || J+M<0 || j1+j2+J<0 || mod(j1+m1,1)~=0 || mod(j2+m2,1)~=0 || mod(J+M,1)~=0
            res=0;
            return
    elseif m1+m2-M~=0 || J<abs(j1-j2) || J>j1+j2
            res=0;
            return
    end



	if J==0
		if j1==j2 && m1==-m2
			res=(-1)^(j1-m1)/sqrt(2*j1+1);
        else
			res=0;
        end
        return
    elseif J==j1+j2 && M==J
		res=1;
        return
    elseif j1==j2 && j1==J/2 && m1==-m2
		res=(factorial(2*j1))^2/(factorial(j1-m1)*factorial(j1+m1)*sqrt(factorial(4*j1)));
        return
    elseif j1==j2 && j1==m1 && j1==-m2
		res=factorial(2*j1)*sqrt((2*J+1)/(factorial(J+2*j1+1)*factorial(2*j1-J)));
        return
    elseif j2==1 && m2==0
		if J==j1+1
			 res=sqrt((j1-M+1)*(j1+M+1)/((2*j1+1)*(j1+1)));
        elseif J==j1
			 res=M/sqrt(j1*(j1+1));
        else
			 res=-sqrt((j1-M)*(j1+M)/((2*j1+1)*j1));
        end
        return
    end

    k=max([0,j1+m2-J,j2-m1-J]):min([j1+j2-J,j1-m1,j2+m2]);

	s=sum((-1).^k./(factorial(k).*factorial(j1+j2-J-k).*factorial(j1-m1-k).*factorial(j2+m2-k).*factorial(J-j2+m1+k).*factorial(J-j1-m2+k))); 

	res=sqrt(((2*J+1)*factorial(J+j1-j2)*factorial(J-j1+j2)*factorial(j1+j2-J))/factorial(j1+j2+J+1))*sqrt(factorial(J+M)*factorial(J-M)*factorial(j1-m1)*factorial(j1+m1)*factorial(j2-m2)*factorial(j2+m2))*s;

end