%This class defines system's dissipation matrix and density matrix using symbolic variables. If numerical values are necessary,
%one can always use function subs(obj.***,symbolic variables,numerical values) and store the result 
%in a new user-defined variable. To find symbolic variables that are being used, use symvar(obj.***).
% 
% The stored properties are:
% -> densityMatrix - the density matrix for which we are solving Bloch equations. To create a dissipation matrix we also need it. 
%
% -> dissipator - the dissipation matrix. Diagonal elements correspond to population transfers due to spontaneous decays. Off-diagonal elements
% are the decoherence terms created by the spontaneous decay only. 
%
% -> branching - a matrix of branching ratios. The initial states are rows, final states are columns, so B(i,f)=b_if means that the branching ratio 
% of the i->f decay is b_if (if state |i> decays with rate G, the decay to f is with rate b_if*G)
%
% -> decayR - a vector storing total decay rates from all the states
%
% 
% The class contains three methods: 
%
% -> contructor Dissipator(n) - to create the object in this class, you have to call this function first. 
% The input is number of states (numerical value). The function creates empty matrices and stores them in appropriate variables. It also creates
% a symbolic n x n density matrix P and saves it. 
%
% -> function addDecay(S_i,S_f,g) - this function adds a decay from initial state |i> to final state |f>. S_i and S_f being indeces are
% numerical values. The variable "g" is the decay rate from |i> to |f> only, so if state |i> decays to several states, and its total decay rate is G,
% and the branching ratio to |f> is b_if, then g=G*b_if. It has to be a symbolic variable.
%
% -> function fromBranching(BR,DR) - if you already have branching ratios and total decay rates of all the states, you can use this function to create the
% dissipator. BR is n x n matrix representing branching ratios in the same form as the property "branching". DR is a vector of al the decay rates, which
% has to be in the same form as property "decayR"
%
% **************************************************************************************************************************************************************
% 
% EXAMPLE - 3-level Lambda system
% 
%                                          |e>
%                                     ------------    State 3
%                                         /  \
%                                        /    \ 
%                                       /      \
%                                  Ga  /        \ Gb
%                                     /         \/
%                                    \/       -------- |b>     State 2
%                  State 1    |a> -------
% 
% Method I - Adding decays one-by-one
% 0. syms Ga Gb real;
% 0. assume(Ga,'positive')
% 0. assume(Gb,'positive')
% 1. d=Dissipator(3);
% 2. d.addDecay(3,1,Ga);
% 3. d.addDecay(3,2,Gb);
% 
% Method II - Using branching ratios
% 0. syms G b31 b32 real;
% 0. assume(G,'positive')
% 1. d=Dissipator(3);
% 2. decays=[0,0,G];
% 3. branchings=[[0,0,0];[0,0,0];[b31,b32,0]]
% 4. d.fromBranching(branchings,decays)
%
% 
% Method I
%                               [0 0 0]
% 1. Creates a 3x3 empty matrix [0 0 0] and saves it to "dissipator" property. Creates another 3x3 empty matrix and saves it to
%                               [0 0 0]                                                                                        [P1_1 P1_2 P1_3]
% to the "branching" property. Creates empty vector [0, 0, 0] and saves it to "decayR". Finally, creates a symbolic 3x3 matrix [P2_1 P2_2 P2_3]
%                                                                                                                              [P3_1 P3_2 P3_3]
% and saves it to "densityMatrix".
% 2. To properly add both diagonal and off-diagonal elements, we create use an operator G=sqrt(Ga)|f><i|, which in matrix representation 
%    [0 0 sqrt(Ga)]
% is [0 0    0    ]. Then, the contributions from the decay are calculated using formula -1/2{G'G,P}+GPG', where P is the density matrix, {}
%    [0 0    0    ]
% represents the anticommutator, and G' is hermitian conjugate (the formula is very general and creates a matrix that in Bloch equations will
%                                                                  [ Ga*P3_3        0      -Ga/2*P1_3]
% perserve the trace). After the operation, the dissipator becomes [    0           0      -Ga/2*P2_3]. At the end, branching ratios and the decay
%                                                                  [-Ga/2*P3_1  -Ga/2*P3_2   -Ga*P3_3]
%                                                   [0 0 0]  
% rates are updated. Branching ratio matrix becomes [0 0 0], and the decay rates vector is [0,0,Ga].
%                                                   [1 0 0]                                                [0 0    0    ]
% 3. Analogical operation is performed. This time the G operator in matrix representation becomes equal to [0 0 sqrt(Gb)], and after the operation
%                                      [     Ga*P3_3            0       -(Ga+Gb)/2*P1_3]                   [0 0    0    ]    [    0          0      0]  
% is performed, the dissipator becomes [        0            Gb*P3_3    -(Ga+Gb)/2*P2_3], while the branching ratios are now [    0          0      0], 
%                                      [-(Ga+Gb)/2*P3_1  -(Ga+Gb)/2*P3_2  -(Ga+Gb)*P3_3]                                     [Ga/(Ga+Gb) Gb/(Ga+Gb) 0]
%
% and the decay rates are [0,0,Ga+Gb].
%                     
% 
% Method II
% 
% 4. When this function is called diagonal elements are first added "by hand". That is, first states indicated by the vector DR have a diagonal element 
% equal to -DR(i)*Pi_i (total decay rate times the population), then, for every initial state, all the final states have a diagonal element equal to
%                                                              [G*b31*P3_3      0          0  ]
% DR(i)*Pi_i*b_if. For the system under consideration we'd get [    0       G*b32*P3_3     0  ]. If the branching ratios are numerical values, the
%                                                              [    0           0      -G*P3_3]
% function checks if the branching ratios have value between 0 and 1, and if for every decaying state they sum to 1 (they have to!). Finally, off-diagonal
%                                                                                                                              [0 0    0   ]
% elements are added. The we can do by simplifying the formula used in Method I. Namely this time we first create a matrix M = [0 0    0   ] and add first
%                                                                                                                              [0 0 sqrt(G)]
% part of the off-diagonal elements (so in more general case the matrix M(i,i)=sqrt(DR(i)) for all i) by calculating -1/2{M^2,P} and adding it to the dissipator.
% The second part is added, by creating matrices Mi such that Mi(j,j)=sqrt(DR(j))*delta(i,j), and then by adding sum over i of Mi*P*Mi to the dissipator.

classdef Dissipator < handle
    properties
        densityMatrix
        dissipator
        branching
        decayR
    end
    
    methods
        function obj=Dissipator(n)
            if isnumeric(n)
                
               P=sym('P',n);
               
               % The symbolic matrix is created in such a way, that the
               % i+1,j+1-th element is piqj
               for i=1:n
                  for j=1:n
                      if verLessThan('matlab', '9.3')
                          P(i,j)=sym(sprintf('p%dq%d(t)',i-1,j-1));
                      else
                          P(i,j)=str2sym(sprintf('p%dq%d(t)',i-1,j-1));
                      end                    
                  end
               end
                
               n=length(P);
               
               D=sym('D',n);
               D(:,:)=sym(0);
               B=sym('B',n);
               B(:,:)=sym(0);
               G=sym('G',[n 1]);
               G(:,:)=sym(0);
               obj.densityMatrix=P;
               obj.dissipator=simplify(D,'Steps',10);
               obj.branching=B;
               obj.decayR=G;
           else
               error('The dimensions has to be numeric value')
           end          
        end

        function obj=addDecay(obj,si,sf,g)
            D=obj.dissipator;
            P=obj.densityMatrix;
            B=obj.branching;
            R=obj.decayR;
            
            if isnumeric(si) && isnumeric(sf)
               if si>length(D) || sf>length(D)
                   error('Indecies cannot exceed size of the hamiltonian')
               end
               if si==sf
                   error('State cannot decay into itself')
               end
               if isnumeric(g)
                   g=sym(g);
               end

            
                n=length(D);
                G=sym('G',n);
                G(:,:)=sym(0);
                G(sf,si)=sqrt(g);
                
                % Adding effects of a single spontaneous decay to the
                % dissipator
                D_new=D-1/2*anticommute(G'*G,P)+G*P*G';

                obj.dissipator=simplify(D_new,'Steps',20);
                
                
                %Update decay rates       
                bold=R(si);
                bnew=R(si)+g;
                R(si)=bnew;
                obj.decayR=R;


                %Update branching ratios
                for i=1:n
                    if i==sf
                        B(si,i)=(B(si,i)*bold+g)/bnew;
                    else
                        B(si,i)=B(si,i)*bold/bnew;
                    end
                end
                obj.branching=B;
           else
               error('Indecies have to be numeric values')
           end
                      
        end
        
        
        function obj=fromBranching(obj,BR,DR)
               D=obj.dissipator;
               P=obj.densityMatrix;
               
               n=length(D);
               
               eps=10^(-3); %Sum of branching ratios must be equal to 1. This determines the tolerance level.
               
               if length(BR)~=n || length(DR)~=n
                   error('Branching ratios and decay rates have to be specified for all the states')
               else
                    if isnumeric(BR)
                      for i=1:n
                        for j=1:n
                            %Check of single branching ratios
                          if BR(i,j)>1 || BR(i,j)<0
                            error('Branching ratios have to be between 0 and 1')
                          end
                        end
                        
                        %Check of sum of branching ratios
                        if (sum(BR(i,:))>1+eps || sum(BR(i,:))<1-eps) && (sum(BR(i,:))>eps || sum(BR(i,:))<-eps)
                          
                          error('Sum of branching ratios for a decay from any state has to be equal to 1 or 0')
                        end
                      end
                    end
                   obj.branching=BR;
                   obj.decayR=DR;
                   
                   %Diagonal elements
                   G=sym('G',n);
                   G(:,:)=sym(0);
                   
                   for i=1:n
                       D(i,i)=-DR(i)*P(i,i);
                       G(i,i)=sqrt(DR(i));
                   end
                   for i=2:n
                       for j=1:i-1
                           D(j,j)=D(j,j)+BR(i,j)*DR(i)*P(i,i);

                       end
                   end   
                   
                   
                   %Off-diagonal
                   D=D-1/2*anticommute(G*G,P);
                   for i=1:n
                       GG=sym('G',n);
                       GG(:,:)=sym(0);
                       GG(i,i)=sqrt(DR(i));
                       D=D+GG*P*GG';
                   end
                   obj.dissipator=simplify(D,'Steps',100);     
               end

        end 
 
    end
end
            
            
            
            
            
function[vec]=to_vector(A)
    k=size(A,1);
    vec=sym('v',[k^2 1]);
    for i=1:k
        vec(k*(i-1)+1:k*i)=A(:,i);
    end
end

function[X]=commute(A,B)
    X=A*B-B*A;
end

function[X]=anticommute(A,B)
    X=A*B+B*A;
end           
       


 
            
            
            
            
            