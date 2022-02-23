%
% This class builds a system of Bloch equations using Hamiltonian and the dissipation matrix. 
%
%
%
% ///////////////////////////////////////////////////////////////////////////////////////
%
% The stored public properties are:
%
% steadyState -> NxN matrix, where (i,j) element corresponds to steady state solution for Pij ((i,j) element of the density matrix). Pii are
% populations.
%
% densityMatrix -> density Matrix P copied from the dissipator object.
%
% evTime -> times at which the Bloch equations where evaluated when solving for the time evolution. It's a vector.
%
% evolution -> NxNxL matrix, where (i,j) element is a vector of length L storing values of Pij at times saved to obj.evTime; i.e. it's the time
% evolution of all elements of the density matrix.
%
% evolutionTr -> NxNxL matrix, where (i,j) element is a vector of length L storing values of Pij in changed basis at times saved to obj.evTime; i.e. it's the time
% evolution of all elements of the density matrix.
%
% eqnsRHS -> N^2 long vector storing the right hand side of Bloch Equations in the vector form, i.e. i[H,P]+D. The order is:
% [(1,1),(1,2),...,(1,N),(2,1)...,(N,N-1),(N,N)]
%        
% equations -> NxN matrix storing full Bloch equations, i.e. dP/dt==i[H,P]+D
%
% equationsVector -> N^2 vector storing full Bloch equations
%
% equationsS -> NxN matrix storing steady state equations, i.e. 0==i[H,P]+D
%
% intTime -> 2-element vector storing the integration time for time evolution. intTime=[initial_time final_time].
%
% initialConditions -> NxN matrix storing initial conditions used for the time evolution.
%
% lastSol -> last solution obtained from numerical integration; used when extending the evolution
%
% optParams -> optimal parameters found through optimization
% 
% optVal -> optimal value found through optimization
%
%
%
%
% /////////////////////////////////////////////////////////////////////////////////////// 
%
% The class contains following methods: 
%
% -> contructor BlochEqns(H,D) - to create the object in this class, you have to call this function first. 
% The input consists of Hamiltonian object H and dissipator object D. The function copies the density matrix from D and sets up all the Bloch
% equations and saves them to appropriate variables. It uses the hamiltonian after the unitary transformation. It uses basic hamiltonian,
% if the transformation wasn't performed.
%
% -> function solveSteady(v) - solves Bloch equations symbolically in the steady state and saves the solution. Note, that the solution usually is feasible for
% small systems.
%
% -> function necessaryVariables() - finds all the symbolic variables in Bloch equations that need to be substituted before solving the time
% evolution. 
%
% -> function evolve(t_initial,t_final,IC,vars) - evolves the Bloch equation using ode15s solver. Evolution happens from 't_initial' to
% 't_final' with initial conditions 'IC'. 'vars' are numerical values of symbolic variables found using 'necessaryVariables()' method. They have to be
% provided in the same order!
%
% -> function plotEvolution() - plots diagonal elements of obj.evolution (populations of states) as function of time
%
% -> function extendEvolution(t_final,Vars) - evolves Bloch equations using the previously obtained solution as the initial condition. Solution is
% concatenated to the already existing one.
%
% -> function changeBasis(vecs) - performs transformation of the density matrix at every found time step (obj.evolution) to a different basis. As argument it takes 'n' 
% (number of states) vectors of length 'n' defining the new basis, combines them in transformation matrix U and performs operation of U'PU. Results
% are stored in 'obj.evolutionTr'.
%
% -> function optimizeParameters(t_ini,t_fin,IC,M,score_pop,varargin) - optimizes parameters using genetic optimization. Can be optimized to obtain high/low populations in
% various states. t_ini and t_fin are the integration time, IC - initial conditions, M - cell matrix with parameters and their domains, e.g.
% M={G,0,1;w_a,-2,2}, score_pop - list with states with regards to which we're optimizing. Finally, there are several options:
% 'Criterion'='minimum' or 'maximum' - minimizing or maximizing 
% 'Popsize'=int>0 - number of random set of parameters at every iteration
% 'Iterations'=int>0 - number of iterations every population is optimized for
% 'Popnumber=int>0 - number of different populations initialized and seperately optimized
% 'Fraction'=1>float>0 - part of population that is carried over to next iteration
% 'Integration'='yes' or 'no' - if 'no' is chosen population in given states at t_fin is optimized; if 'yes' is chosen, integral of the
% population as a function of time from t_ini to t_fin is optimized
% 'Switching'='yes' or 'no' - if 'yes' is chosen every
% t_step=(t_fin-t_ini)/no_switches set of parameters used for optimization is being turned off/on
% 'NoSwitches'=int>0 - how many times parameters are switched in the integration time
% 'SwitchingParameters'=matrix(a,b) - parameters showing when to switch on/off optimizable parameters; matrix should have at most 'NoSwitches'
% number of rows and number of columns equal to number of optimized parameters; in a specific row '0' indicates parameter being turned off,
% any other number indicates parameter turned on
% Default={'Criterion':'maximum'; 'Popsize':20; 'Iterations':10; 'Popnumber':3; 'Fraction':0.3; 'Integration': 'no'; 'Switching':'no';
% 'NoSwitches':0; 'SwitchingParameters':[1,1,...,1] (length=number of optimized parameters)}
%
%
%
% /////////////////////////////////////////////////////////////////////////////////////// 




classdef BlochEqns < matlab.mixin.Copyable%handle
    
    
    properties
        steadyState
        densityMatrix
        evTime
        evolution
        evolutionTr
        eqnsRHS
        equations
        equationsVector
        equationsS
        intTime
        lastSol
        initialConditions
        optParams
        optVal
    end
    
    
    methods
        function obj=BlochEqns(H,D)
            import bloch.Hamiltonian
            import bloch.Dissipator

            if ~isa(H,'Hamiltonian') || ~isa(D,'Dissipator')
                error('You have to provide hamiltonian and dissipators as objects of their own class')
            else
                n=length(H.transformed);
                if n==0
                    n=length(H.hamiltonian);
                end
                if length(D.dissipator)~=n
                    error('Dissipator and hamiltonian have to be of the same size')
                else
                    P=D.densityMatrix;
                    
                    t=sym('t','real');
                    
                    
                    Ham=H.transformed;
                    if isempty(Ham)
                        Ham=H.hamiltonian;
                    end
                    L=D.dissipator;
                    
                    Final_matrix=simplify(-1i*commute(Ham,P)+L,'Steps',50); %Master equation
                    Eqns=diff(P,t)==Final_matrix;
                    Eqns0= zeros(n)==Final_matrix;
                    Eqns0=to_vector(Eqns0);
                    Eqns0=[Eqns0;1==trace(P)]; %Steady state equations need trace(P)==1 added to be fully solved
                    
                    obj.equations=Eqns;
                    obj.equationsVector=to_vector(Eqns);
                    obj.equationsS=Eqns0;

                    obj.eqnsRHS=to_vector(Final_matrix);
                    obj.densityMatrix=P;
                    
                    disp('Remember to call "necessaryVariables" method before solving the time evolution of the system!')
                end
            end
            
        end
        
        
        function obj=solveSteady(obj)
            
            n=length(obj.densityMatrix);
            Eqs=obj.equationsS ;%Solver requires equations to be in a vector
            
            Variables=to_vector(obj.densityMatrix); 
            
            Y=sym('Y', [n^2 1]);

            Eqs=subs(Eqs,Variables,Y); %They also have to be solved by substituting density matrix terms with dummy variables
            
            sol=solve(Eqs,Y);
            
            %We put solutions back into a matrix form
            SM=sym('SM',n);
            for i=1:n
                for j=1:n
                    SM(i,j)=simplify(sol.(['Y' num2str((i-1)*n+j)]),'Steps',20);
                end
            end
            obj.steadyState=simplify(SM,'Steps',100);
        
        end
        
        
        function v=necessaryVariables(obj)
            v=symvar(obj.eqnsRHS);
            syms t;
            disp('To solve the time evolution of the system you have to provide numerical values for the following variables:')
            disp(setdiff(v,t))
            disp('The order has to be the same')
            v=setdiff(v,t);
        end
        
        
        function obj=evolve(obj,t_initial,t_final,IC,Vars)
            if isnumeric(t_final)
                if t_final<=t_initial
                    error('Time of the evolution must be greater than 0')
                else
                    obj.intTime=[t_initial,t_final];
                    n=length(obj.densityMatrix);
                    if length(IC)~=n
                        error('Please provide initial conditions for all states.')
                    else
                        %All the equations and initial conditions have to
                        %be vectorized
                        obj.initialConditions=IC;
                        IC=to_vector(IC);
                      
                        Variables=to_vector(obj.densityMatrix);  
                        
                        Y=sym('Y', [n^2 1]); %Dummy variables
                        t=sym('t','real');

                        Eqs=obj.eqnsRHS;                    
                        v=symvar(Eqs);
                        Add_vars=setdiff(v,t);
                        
                        
                        F=subs(Eqs,[transpose(Variables),Add_vars],[transpose(Y),Vars]);

                        FF=matlabFunction(F,'vars',{t,Y}); %Matlab requires us to create a matlabFunction object from our equations, in order to solve ODEs
                        
                        %AbsTol and RelTol are tolerances that determine
                        %how precise our solution is. MaxStep given here
                        %limits the solver to take minimum of 1000 steps.
                        option=odeset('AbsTol',1e-6,'RelTol',1e-5,'MaxStep',(t_final-t_initial)/1000);
      
                        solution=ode45(FF,[t_initial t_final],double(IC),option);
                                               
                        obj.lastSol=solution;
                     
                        obj.evTime=solution.x(1,:); %Times at which solution was found
                        
                        %Solutions are put into a matrix
                        EV=zeros(n,n,length(solution.x(1,:))); 
                        for i=1:n
                            for j=1:n
                                EV(j,i,:)=solution.y((i-1)*n+j,:);
                            end
                        end

                        obj.evolution=EV;
                    end
                end             
            else
                error('Time of the evolution must be numeric value')
            end    
        end
        
        
        function x=getForce(obj,H,time)
           import bloch.Hamiltonian
           if isempty(obj.evolution)
               error('You need to obtain time evolution of the density matrix first')
           end
           
           if ~isa(H,'Hamiltonian')
                error('You have to provide hamiltonian as object of its own class')
           end
           
           if time<obj.intTime(1) || time >obj.intTime(2)
               error('Time at which force is evaluated should be withing limits of the integration time used')
           end
           
           syms t real;
           
           Hvars=symvar(H.hamGradient);
           
           if length(Hvars)>2
               error('Gradient should have time [t] as its only symbolic variable')
           elseif ~ismember(t,Hvars)
               
               error('Gradient should have time [t] as its only symbolic variable')
           end
           
           mn=obj.intTime(2)-obj.intTime(1);
           ind=1;
           for i=1:length(obj.evTime)
              if abs(obj.evTime(i)-time)<mn
                  mn=abs(obj.evTime(i)-time);
                  ind=i;
              else
                  break
              end
           end
           
           gradH=subs(H.hamGradient,t,obj.evTime(ind));
           
           x=-double(real(vpa(simplify(trace(squeeze(obj.evolution(:,:,ind))*gradH)))));
            
        end
        
        function x=averageForce(obj,H,time_ini,time_fin,step)
           import bloch.Hamiltonian
           if isempty(obj.evolution)
               error('You need to obtain time evolution of the density matrix first')
           end
           
           if ~isa(H,'Hamiltonian')
                error('You have to provide hamiltonian as object of its own class')
           end
           
           if time_ini<obj.intTime(1) || time_ini>obj.intTime(2)
               error('Time at which force is evaluated should be withing limits of the integration time used')
           end
           
           if time_ini==time_fin
               error('Provided times have to be different')
           end
           
           if time_fin<obj.intTime(1) || time_fin>obj.intTime(2)
               error('Time at which force is evaluated should be withing limits of the integration time used')
           end
           
           if time_fin<time_ini
               temp=time_fin;
               time_fin=time_ini;
               time_ini=temp;
           end
           
           
           syms t real;
           
           Hvars=symvar(H.hamGradient);
           
           if length(Hvars)>2
               error('Gradient should have time [t] as its only symbolic variable')
           elseif ~ismember(t,Hvars)
               error('Gradient should have time [t] as its only symbolic variable')
           end
           
           mn=obj.intTime(2)-obj.intTime(1);
           ind_i=1;           
           for i=1:length(obj.evTime)
              if abs(obj.evTime(i)-time_ini)<mn
                  mn=abs(obj.evTime(i)-time_ini);
                  ind_i=i;
              else
                  break
              end
           end
           
           
           
           mn=obj.intTime(2)-obj.intTime(1);
           ind_f=length(obj.evTime);           
           for i=length(obj.evTime):-1:1
              if abs(obj.evTime(i)-time_fin)<mn
                  mn=abs(obj.evTime(i)-time_fin);
                  ind_f=i;
              else
                  break
              end
           end
           
           Force=zeros([round((ind_f-ind_i)/step)+1,1]);
           k=1;
           for i=ind_i:step:ind_f
              gradH=subs(H.hamGradient,t,obj.evTime(i));          
              x=-real(trace(squeeze(obj.evolution(:,:,i))*gradH));
              Force(k)=x;
              k=k+1;
           end
           x=double(trapz(obj.evTime(ind_i:step:ind_f),Force)/(time_fin-time_ini));           
           
        end
        
        function obj=changeBasis(obj,varargin)
           Eqs=obj.evolution;
           ne=length(obj.densityMatrix);
           n=nargin-1;

           if n~=ne
               error('Number of vectors provided has to be equal to the number of states')
           end
           
           U=[];
           
           for i=1:n
               vec=varargin{i};

               if length(vec)~=ne
                
                   error('All vectors have to have the same length as number of states')
                   
               end
               s=size(vec);
               
               if (s(1)==1 && s(2)==1 && ne>1) || length(s)>2 || (s(1)>1 && s(2)>1)
                   error('You have to provide vectors.')
               end
               if s(1)==1
                   vec=vec.';
               end
               
               U=[U,vec];
               
           end
           
           EqsTr=zeros(size(Eqs));
           
           for i=1:size(Eqs,3)
               EqsTr(:,:,i)=U'*squeeze(Eqs(:,:,i))*U; %Basis is changed at every time step
           end
           
           obj.evolutionTr=EqsTr;
                       
        end
       
        
        function extendEvolution(obj,t_final,Vars)
            if isnumeric(t_final)
                t_initial=obj.intTime(2);
                    if t_final<=t_initial
                        error('Time of the evolution must be greater than 0')
                    else
                        
                        n=length(obj.densityMatrix);

                            Variables=to_vector(obj.densityMatrix);  

                            Y=sym('Y', [n^2 1]);
                            t=sym('t','real');

                            Eqs=obj.eqnsRHS;                    
                            v=symvar(Eqs);
                            Add_vars=setdiff(v,t);

                            F=subs(Eqs,[transpose(Variables),Add_vars],[transpose(Y),Vars]);

                            FF=matlabFunction(F,'vars',{t,Y});

                            solution=odextend(obj.lastSol,FF,t_final);
                            
                            obj.lastSol=solution;

                            obj.evTime=solution.x(1,:);
                            
                            EV=zeros(n,n,length(solution.x(1,:)));
                            for i=1:n
                                for j=1:n
                                    EV(j,i,:)=solution.y((i-1)*n+j,:);
                                end
                            end

                            obj.evolution=EV;
                    end             
                else
                    error('Time of the evolution must be numeric value')
            end
        end
        
        
        function plotEvolution(obj)
            figure
            for i=1:length(obj.densityMatrix)
                plot(obj.evTime(1,:),squeeze(obj.evolution(i,i,:)))
                hold on
                legendInfo{i}=['State ' num2str(i)];  %Generic names are provided
            end
            xlabel('Time')
            ylabel('Population')
            xlim([obj.evTime(1,1),obj.evTime(1,end)])
            legend(legendInfo)
            drawnow
        end

        %Genetic optimization function. The parameters to be optimized have
        %to be provided as a cell matrix
        %'M'={'Var1',lower_bound,upper_bound;'Var2',...}. 'score_pop'
        %should be a list of indices. The variables to be optimized given
        %in 'M' have to be the same as ones provided by
        %'necessaryVariables' method, and have to be provided in the same
        %order.
        function obj=optimizeParameters(obj,t_ini,t_fin,IC,M,score_pop,varargin)

            rng shuffle;

            sz=size(M);
            if sz(2)~=3
                error('You have to provide lower and upper bounds for all parameteres being optimized. Please provide all the parameters and bounds in form of Nx3 matrix.');
            end
            if ~iscell(M)
                error('The matrix with parameters has to be defined as a cell');
            end

            N=sz(1);
            Domain=zeros(N,2); %Bounds for optimized parameters. 
            
            %Checking if parameters for optimization are viable
            nec_vars=obj.necessaryVariables();
            for i=1:N
                if isnumeric(M{i,1})
                    error('The first column has to contain symbolic parameters that are being optimized.');
                end
                if ~ismember(M{i,1},nec_vars)
                    error('You have to provide parameters that have not been assigned numerical values yet');
                end
                if M{i,1}~=nec_vars(i)
                    error('Parameters have to be provided in the same order as given by the "necessaryVariables()" function');
                end

                if ~isnumeric(M{i,2}) || ~isnumeric(M{i,3})
                    error('Lower or upper bounds have to be numerical values');
                end
                if M{i,2}>=M{i,3}
                    error('Lower bound has to be smaller than the upper bound');
                end

                Domain(i,1)=M{i,2};
                Domain(i,2)=M{i,3};
            end
            
            %Indices of states to be optimized
            if length(score_pop)<1
                error('You have to provide at least one population state with respect to which optimize');
            elseif length(score_pop)>=length(obj.densityMatrix)
                error('You can optimize with respect to at most n-1 populations');
            end
            
            %Validity of indices
            for k=score_pop
                if ~isnumeric(k)
                    error('Indecies have to be numeric');
                elseif k>length(obj.densityMatrix) || k<1
                    error('Indecies cannot bigger than number of states or smaller than 1');
                elseif k-ceil(k)~=0
                    error('Indecies have to be integers');
                end
            end

            %Validating functions for the options
            ScalarValidity=@(x) isnumeric(x) && isscalar(x) && (x>0);
            
            ScalarValidityF=@(x) isnumeric(x) && isscalar(x) && (x>0) && (x<1);
            
            MatrixValidity=@(x) isnumeric(x) && ismatrix(x) && validate_matrix(x,N);
            
            IntegerValidity=@(x) isnumeric(x) && isscalar(x) && (x>0) && mod(x,1)==0;

            expectedCriteria={'minimum','maximum'};
            StringValidityC=@(x) any(validatestring(x,expectedCriteria));
            
            expectedInt={'yes','no'};
            StringValidityI=@(x) any(validatestring(x,expectedInt));

            p=inputParser;

            %Default controls
            default_crit='maximum'; %Criterion for optimization
            default_popsize=20;  %Single population size (number of parameter vectors)
            default_elitefrac=0.3; %Fraction of population that is allowed to have offspring
            default_iterations=10; %Optimization iterations 
            default_populations=3; %Number of times new random populations are initialized nad optimization repeated
            default_int='no'; %Whether or not we optimize absolute value at t_final or integral from t_initial to t_final
            default_sw='no'; %Whether or not we should switch some parameters on and off when solving for time evolution
            default_nosw=0; %Number of times switching happens
            default_swpar=ones(1,N);
            
            %Associating validation functions with specific options
            addParameter(p,'Criterion',default_crit,StringValidityC);
            addParameter(p,'Popsize',default_popsize,IntegerValidity);
            addParameter(p,'Fraction',default_elitefrac,ScalarValidityF);
            addParameter(p,'Iterations',default_iterations,IntegerValidity);
            addParameter(p,'Popnumber',default_populations,IntegerValidity);
            addParameter(p,'Integration',default_int,StringValidityI);
            addParameter(p,'Switching',default_sw,StringValidityI);
            addParameter(p,'NoSwitches',default_nosw,IntegerValidity);
            addParameter(p,'SwitchingParameters',default_swpar,MatrixValidity);
    
            parse(p,varargin{:}); %option parser

            criterion=p.Results.Criterion;
            popsize=p.Results.Popsize;
            elite_frac=p.Results.Fraction;
            iterations=p.Results.Iterations;
            distinct_pops=p.Results.Popnumber;
            integration=p.Results.Integration;
            switching=p.Results.Switching;
            no_switches=p.Results.NoSwitches;
            switching_parameters=p.Results.SwitchingParameters;
            
            
            diff_switches=length(switching_parameters(:,1));
            
            if strcmp(switching,'yes')
                if no_switches<=1
                    error('When optimizing evolution with parameter switching, please provide the number of switches to consider (option "NoSwitches")');
                end
                if diff_switches>no_switches || diff_switches<1 
                    error('Information which prameter is to be switched on or off has to be provided for all switches. Therefore, if swtiching is enabled, please provide a at least two vectors of 0s and 1s, and not more than number of switches');
                end
            end
            
            
            elite=int8(elite_frac*popsize); %Absolute number of examples allowed to have offspring and move to the next iteration
            final_elite=max(int8(popsize/2/distinct_pops),3); %After multiple rounds (specified in options) of optimizations are performed, 
            %top 'final_elite' number of examples from the last iteration is copied to the final population and allowed to have new offspring


            % Probabilities of mutation and crossbreeding
            probs=zeros(N+1,1);

            probs(1)=0.2;
            if N>4
                val=0.15/(N-4);
                probs(N+1)=0.6;
                for k=2:N
                    if k==2
                        probs(k)=0.35;
                    elseif k>2 && k<5
                        probs(k)=probs(k-1)+0.05;
                    else
                        probs(k)=probs(k-1)+val;
                    end
                end
            elseif N==3
                probs(2)=0.35;
                probs(3)=0.4;
                probs(4)=0.41;
            elseif N==2
                probs(2)=0.35;
                probs(3)=0.36;
            elseif N==4
                probs(2)=0.35;
                probs(3)=0.4;
                probs(4)=0.45;
                probs(5)=0.46;
            end


            best=zeros(1,N+1);
            if strcmp(criterion,'minimum')
                 best(1,N+1)=10^6;
            end

            %Initialization

            for l=1:distinct_pops
                fprintf('Population %d \n',l);
                pop=zeros(popsize,N,'single');

                for i=1:popsize
                    for j=1:N
                        pop(i,j)=rand*(Domain(j,2)-Domain(j,1))+Domain(j,1);
                    end
                end


                %Optimization
                for i=1:iterations

                   fprintf('Iteration %d \n',i);

                   scores=zeros(popsize,N+1,'single');
                   scores(:,1:N)=pop;
                   for j=1:popsize
                         vars=zeros(1,N);
                         for k=1:N
                             vars(k)=pop(j,k);
                         end
                         
                         
                         if strcmp(switching,'no')
                            obj.evolve(t_ini,t_fin,IC,vars);
                         else
                             t_step=(t_fin-t_ini)/no_switches;
                             for sw=1:no_switches
                                
                                 row=mod(sw,diff_switches);
                                 if row==0
                                     row=diff_switches;
                                 end
                                 var_sw=vars;
                                 
                                 for i=1:N
                                     if switching_parameters(row,i)==0
                                         var_sw(i)=0;
                                     end
                                 end
                                 
               
                                if sw==1
                                    obj.evolve(t_ini,t_step,IC,var_sw); %We solve equations for every vector of parameters in the population
                                else
                                    obj.intTime=[(sw-2)*t_step,(sw-1)*t_step];
                                    obj.lastSol=prevSol;
                                    obj.extendEvolution(t_step*sw,var_sw);
                                end

                                prevSol=obj.lastSol;
                             end
                         end
                         
                         %Parameter vectores are ordered by the score the
                         %got, which here is simply min/max of the sum of
                         %population at provided indices at t_final (or sum
                         %of integrals if chosen so in the options)
                         scores(j,N+1)=0;
                         if strcmp(integration,'yes')
                             Ys=zeros(1,length(obj.evTime(:)));
                             for ind=score_pop
                                 Ys(1,:)=Ys(1,:)+single(real(squeeze(obj.evolution(ind,ind,:))'));
                             end
                             scores(j,N+1)=trapz(obj.evTime(:),Ys(:));
                         else
                             for ind=score_pop
                                 scores(j,N+1)=scores(j,N+1)+single(real(obj.evolution(ind,ind,end)));
                             end
                         end

                         obj.evolution(:,:,:)=[];
                         obj.evTime(:)=[];
                   end

                   if strcmp(criterion,'maximum')
                       opt='descend';
                   elseif strcmp(criterion,'minimum')
                       opt='ascend';
                   end

                   [Y,I]=sort(scores(:,N+1),opt);
                   clear Y;
                   scores=scores(I,:);
                   clear I;

                   ranked=scores(:,1:N);
                   
                   if strcmp(criterion,'maximum')
                       if scores(1,N+1)>best(1,N+1)
                           best(1,:)=scores(1,:);
                           disp(best); %best parameters and the score are shown
                       end
                   elseif strcmp(criterion,'minimum')
                       if scores(1,N+1)<best(1,N+1)
                           best(1,:)=scores(1,:);
                           disp(best); %best parameters and the score are shown
                       end
                   end
                   

                   clear scores;

                   pop=ranked(1:elite,:);


                   %Loop that uses 'elite' number of top parameter vectors
                   %to replenish the population using either: mutation or
                   %crossbreeding. Mutation can happen at a single, few or
                   %all sites. 
                   while size(pop(:,1))<popsize
                      pr=rand;
    
                      if pr<probs(N) %Mutation
                          if N==1 || pr<probs(1)
                              c1=randi(N);
                              C=[c1];
                          else
                            for k=2:N              
                                if pr<probs(k) && pr>probs(k-1)
                                    if k<N/2+1
                                      C=randi(N,1,k);
                                      C=unique(C);
                                      while length(C)<k
                                          ck=randi(N);
                                          if ~ismember(ck,C)
                                              C=[C,ck];
                                          end
                                      end
                                    elseif k==N
                                        C=1:N;
                                    else
                                        C=1:N;
                                        D=randi(N,1,N-k);
                                        D=unique(D);
                                        while length(D)<N-k
                                            ck=randi(N);
                                            if ~ismember(ck,D)
                                                D=[D,ck];
                                            end
                                        end
                                        C(D)=[];
                                    end
                                   
                                  if length(C)==k
                                      break
                                  end
                                end
                            end
                          end

                          vec=pop(randi(elite),:);

                          for c=C                                           
                               vec(c)=vec(c)+(rand*2-1)*0.1*(Domain(c,2)-Domain(c,1));
                               if vec(c)<Domain(c,1)
                                   vec(c)=Domain(c,1);
                               elseif vec(c)>Domain(c,2)
                                   vec(c)=Domain(c,2);
                               end
                          end


                          pop=[pop;vec];

                          clear vec;

                      elseif k<probs(N+1) %Complete randomization. Initializes a new parameter vector.
                          vec=zeros(1,N);
                          for ind=1:N
                            vec(ind)=rand*(Domain(ind,2)-Domain(ind,1))+Domain(ind,1);
                          end

                          pop=[pop;vec];
                          clear vec;    
                      else %Crossbreeding
                          c1=randi(elite);
                          c2=randi(elite);
                          vec1=pop(c1,:);
                          vec2=pop(c2,:);
                          r=randi(N);
                          vec(1:r)=vec1(1:r);
                          vec(r+1:N)=vec2(r+1:N);

                          pop=[pop;vec];

                          clear vec vec1 vec2;
                      end

                   end


                end

                final_population((l-1)*final_elite+1:l*final_elite,:)=ranked(1:final_elite,:); %After all iterations, top parameter vectors are added

                clear pop ranked;
            end
            
            
            
            %Using 'final_population', the population for last optimization
            %run is replenished.
            if distinct_pops>1
                pop=final_population;
                clear final_population;
                
               
                if size(pop(:,1))>popsize
                    pop=pop(randperm(size(pop(:,1))),:);
                    pop=pop(1:popsize,:);
                end

                while size(pop(:,1))<popsize
                          pr=rand;

                          if pr<probs(N) 
                              if N==1 || pr<probs(1)
                                  c1=randi(N);
                                  C=[c1];
                              else
                                for k=2:N              
                                    if pr<probs(k) && pr>probs(k-1)
                                        if k<N/2
                                          C=randi(N,1,k);
                                          C=unique(C);
                                          while length(C)<k
                                              ck=randi(N);
                                              if ~ismember(ck,C)
                                                  C=[C,ck];
                                              end
                                          end
                                        elseif k==N
                                            C=1:N;
                                        else
                                            C=1:N;
                                            D=randi(N,1,N-k);
                                            D=unique(D);
                                            while length(D)<N-k
                                                ck=randi(N);
                                                if ~ismember(ck,D)
                                                    D=[D,ck];
                                                end
                                            end
                                            C(D)=[];
                                        end
                                      if length(C)==k
                                          break
                                      end
                                    end
                                end
                              end

                              vec=pop(randi(elite),:);

                              for c=C                                           
                                   vec(c)=vec(c)+(rand*2-1)*0.1*(Domain(c,2)-Domain(c,1));
                                   if vec(c)<Domain(c,1)
                                       vec(c)=Domain(c,1);
                                   elseif vec(c)>Domain(c,2)
                                       vec(c)=Domain(c,2);
                                   end
                              end


                              pop=[pop;vec];

                              clear vec;

                          elseif k<probs(N+1)
                              vec=zeros(1,N);
                              for ind=1:N
                                vec(ind)=rand*(Domain(ind,2)-Domain(ind,1))+Domain(ind,1);
                              end

                              pop=[pop;vec];
                              clear vec;    
                          else
                              c1=randi(elite);
                              c2=randi(elite);
                              vec1=pop(c1,:);
                              vec2=pop(c2,:);
                              r=randi(N);
                              vec(1:r)=vec1(1:r);
                              vec(r+1:N)=vec2(r+1:N);

                              pop=[pop;vec];

                              clear vec vec1 vec2;
                          end

                end

                
                
                
                %Final optimization using winners of previous optimization
                %rounds.
                 disp('Optimization - final population')

                 for i=1:iterations

                       fprintf('Iteration %d \n',i);

                       scores=zeros(popsize,N+1,'single');
                       scores(:,1:N)=pop;
                       for j=1:popsize
                             vars=zeros(1,N);
                             for k=1:N
                                 vars(k)=pop(j,k);
                             end

                             if strcmp(switching,'no')
                                obj.evolve(t_ini,t_fin,IC,vars);
                             else
                                 t_step=(t_fin-t_ini)/no_switches;
                                 for sw=1:no_switches

                                     row=mod(sw,diff_switches);
                                     if row==0
                                         row=diff_switches;
                                     end
                                     var_sw=vars;

                                     for i=1:N
                                         if switching_parameters(row,i)==0
                                             var_sw(i)=0;
                                         end
                                     end

                                     

                                    if sw==1
                                        obj.evolve(t_ini,t_step,IC,var_sw);
                                    else
                                        obj.intTime=[(sw-2)*t_step,(sw-1)*t_step];
                                        obj.lastSol=prevSol;
                                        obj.extendEvolution(t_step*sw,var_sw);
                                    end

                                    prevSol=obj.lastSol;
                                 end
                             end

                             scores(j,N+1)=0;

                             for ind=score_pop
                                 scores(j,N+1)=scores(j,N+1)+single(real(obj.evolution(ind,ind,end)));
                             end

                             obj.evolution(:,:,:)=[];
                             obj.evTime(:)=[];
                       end

                       if strcmp(criterion,'maximum')
                           opt='descend';
                       elseif strcmp(criterion,'minimum')
                           opt='ascend';
                       end

                       [Y,I]=sort(scores(:,N+1),opt);
                       clear Y;
                       scores=scores(I,:);
                       clear I;

                       ranked=scores(:,1:N);

                       if strcmp(criterion,'maximum')
                           if scores(1,N+1)>best(1,N+1)
                               best(1,:)=scores(1,:);
                               disp(best);
                           end
                       elseif strcmp(criterion,'minimum')
                           if scores(1,N+1)<best(1,N+1)
                               best(1,:)=scores(1,:);
                               disp(best);
                           end
                       end

                       clear scores;

                       pop=ranked(1:elite,:);



                       while size(pop(:,1))<popsize
                          pr=rand;

                          if pr<probs(N)
                              if N==1 || pr<probs(1)
                                  c1=randi(N);
                                  C=[c1];
                              else
                                for k=2:N              
                                    if pr<probs(k) && pr>probs(k-1)
                                        if k<N/2+1
                                          C=randi(N,1,k);
                                          C=unique(C);
                                          while length(C)<k
                                              ck=randi(N);
                                              if ~ismember(ck,C)
                                                  C=[C,ck];
                                              end
                                          end
                                        elseif k==N
                                            C=1:N;
                                        else
                                            C=1:N;
                                            D=randi(N,1,N-k);
                                            D=unique(D);
                                            while length(D)<N-k
                                                ck=randi(N);
                                                if ~ismember(ck,D)
                                                    D=[D,ck];
                                                end
                                            end
                                            C(D)=[];
                                        end

                                      if length(C)==k
                                          break
                                      end
                                    end
                                end
                              end

                              vec=pop(randi(elite),:);

                              for c=C                                           
                                   vec(c)=vec(c)+(rand*2-1)*0.1*(Domain(c,2)-Domain(c,1));
                                   if vec(c)<Domain(c,1)
                                       vec(c)=Domain(c,1);
                                   elseif vec(c)>Domain(c,2)
                                       vec(c)=Domain(c,2);
                                   end
                              end


                              pop=[pop;vec];

                              clear vec;

                          elseif k<probs(N+1)
                              vec=zeros(1,N);
                              for ind=1:N
                                vec(ind)=rand*(Domain(ind,2)-Domain(ind,1))+Domain(ind,1);
                              end

                              pop=[pop;vec];
                              clear vec;    
                          else
                              c1=randi(elite);
                              c2=randi(elite);
                              vec1=pop(c1,:);
                              vec2=pop(c2,:);
                              r=randi(N);
                              vec(1:r)=vec1(1:r);
                              vec(r+1:N)=vec2(r+1:N);

                              pop=[pop;vec];

                              clear vec vec1 vec2;
                          end

                       end


                 end
            end

             obj.optParams=best(1,1:N);
             obj.optVal=best(1,N+1);

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


function[bool]=validate_matrix(X,n)
    bool=1;
    sz=size(X);
    if sz(2)~=n
        bool=0;
        return
    end
    X=X(:);
    for i=1:length(X)
        if X(i)~=0 && X(i)~=1
            bool=0;
            return
        end
    end

end




    
