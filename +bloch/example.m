clear
import bloch.*

warning('off','all');


%-----------------------------HAMILTONIAN---------------------------------%

%Define variables
syms wa wb we Wae wae Wbe wbe dae dbe real; %(Rabi rate W in general can be complex)

%Create the hamiltonian
h=Hamiltonian(3);

%Define states of your energies
h.addEnergies([wa,wb,we]);

%Add appropriate coupling for your system (Rabi rates can in principle be time-dependent functions)
h.addCoupling(1,3,Wae,wae);
h.addCoupling(2,3,Wbe,wbe);

%Define detunings for added transitions 
h.defineStateDetuning(1,3,dae);
h.defineStateDetuning(2,3,dbe);

%Define zero energy level
h.defineZero(we);

%Move to a different frame and eliminate time dependencies
h.unitaryTransformation();


disp('Hamiltonian before and after unitary transformation')
disp(h.hamiltonian)
disp(h.transformed)

disp('Matrix used for unitary transformation')
disp(h.transMatrix)



%-----------------------------DISSIPATOR----------------------------------%

%Define variables
syms Ga Gb real;

%One can make assumptions about them if necessary
assume(Ga,'positive')
assume(Gb,'positive')

%Define the Dissipator
d=Dissipator(3);

%Add existing decay paths
d.addDecay(3,1,Ga);
d.addDecay(3,2,Gb);


disp('Dissipator')
disp(d.dissipator)
disp('Decay Rates')
disp(d.decayR)
disp('Branching Ratios')
disp(d.branching)


%If helpful, graph the system
h.createGraph(d);
h.plotGraph();



%----------------------------MASTER EQUATION------------------------------%

%Create the equations
eq=BlochEqns(h,d);

disp('Optical Bloch Equations')
disp(eq.equationsVector)



%----------------------------STEADY STATE---------------------------------%

%Steady state solution, if needed
eq.solveSteady();

disp('Steady state solution for the excited state population')
pee=eq.steadyState(3,3);
disp(pee)



%----------------------------TIME EVOLUTION-------------------------------%

%Before solving ODEs, display variables needed to be substituted
eq.necessaryVariables();

%Define initial conditions
IC=zeros(3);  
IC(1,1)=1;

%Define numerical values for necessary substitutions
gamma_a=0.6;
gamma_b=0.4;
Rabi_a=1;
Rabi_b=1;
delta_a=0;
delta_b=0;
t_initial=0;
t_final=20;

%One can substitute frequently used variables beforehand
eq.eqnsRHS=subs(eq.eqnsRHS,[Ga,Gb],[gamma_a,gamma_b]);

%Solve the master equation
eq.evolve(t_initial,t_final,IC,[Rabi_a,Rabi_b,delta_a,delta_b]);

%And plot it
eq.plotEvolution();

% Extend the evolution, if required
t_final_new=100;
eq.extendEvolution(t_final_new,[Rabi_a,10,delta_a,delta_b])
eq.plotEvolution()




%----------------------------OPTIMIZATION---------------------------------%

%Define domains (upper and lower bounds) for the optimized parameters
Params={Wae,0,1;Wbe,0,1;dae,-5,5;dbe,-5,5};

% Find highest possible population in the excited state
eq.optimizeParameters(t_initial,t_final,IC,Params,[3],'Integration','no','Popsize',10,'Iterations',10,'Popnumber',1);

%Look up optimum parameters
fprintf('Wae = %.2f \n',eq.optParams(1));
fprintf('Wbe = %.2f \n',eq.optParams(2));
fprintf('dae = %.2f \n',eq.optParams(3));
fprintf('dbe = %.2f \n',eq.optParams(4));

% Evolution given optimal parameters
eq.evolve(t_initial,t_final,IC,eq.optParams);
eq.plotEvolution();





