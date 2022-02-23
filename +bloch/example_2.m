clear
import bloch.*

warning('off','all');

%A more complicated system - 4 levels with dark and bright states.

%-----------------------------HAMILTONIAN---------------------------------%

%Define variables and assumptions
syms wg w0 we D real;
syms W_L w_L W_m w_m  real;
syms d_L d_m real;
syms a b A B real;
assume(W_L>0)
assume(W_m>0)

%Create the hamiltonian
H=Hamiltonian(4);

%Define states of your energies
H.addEnergies([w0,wg,wg+D,we]);

%Add appropriate coupling for your system (Rabi rates can in principle be time-dependent functions)
H.addCoupling(1,2,W_m*a,w_m);
H.addCoupling(1,3,W_m*b,w_m);
H.addCoupling(2,4,W_L*A,w_L);
H.addCoupling(3,4,W_L*B,w_L);

%Define zero energy level
H.defineZero(we);

%Move to a different frame and eliminate time dependencies
H.unitaryTransformation();

%Define detunings for added transitions 
H.defineEnergyDetuning(w0,wg,d_m,w_m);
H.defineEnergyDetuning(wg,we,d_L,w_L);


disp('Hamiltonian before and after unitary transformation')
disp(H.hamiltonian)
disp(H.transformed)


%-----------------------HAMILTONIAN MANIPULATION--------------------------%

%Find dark states
H.findDarkStates([2,3],[4]); %2,3 are ground states, 4 is decaying excited state, while 1 is an auxiliary state,
                             %which can be treated as an excited state without decay. We are finding dark states
                             %with respect to transition to the excited state 4.

%Display them
disp('Dark state')
disp(H.darkStates)

%And display the bright states
disp('Bright state')
disp(H.brightStates)

%Move to a basis of dark and bright states. First, create vectors:
u1=[1,0,0,0];
u2=1/sqrt(A^2+B^2)*[0,A,B,0]; %Bright
u3=1/sqrt(A^2+B^2)*[0,B,-A,0]; %Dark
u4=[0,0,0,1];

%Change basis and store the hamiltonian in the original basis beforehand
Htr=H.transformed;
H.changeBasis(u1,u2,u3,u4);

%Display the new hamiltonian (compare H(2,4) and H(3,4))
disp('Hamiltonian in the new basis')
disp(H.transformed)


%If detunings for transition to the auxiliary states are big, we can
%eliminate it using adiabatic elimination:
H.adiabaticElimination([2,3],[1]);

%Display resulting effective hamiltonian
disp('Hamiltonian after adiabatic elimination')
disp(H.adiabaticHamiltonian)


%-----------------------------DISSIPATOR----------------------------------%

%Define variables
syms g real;

%One can make assumptions about them if necessary
assume(g,'positive')


%Define the Dissipator
L=Dissipator(4);

%Add existing decay paths
L.addDecay(4,3,g/2);
L.addDecay(4,2,g/2);


disp('Dissipator')
disp(L.dissipator)


%----------------------------MASTER EQUATION------------------------------%

%Let's bring back the original basis in the hamiltonian
H.transformed=Htr;

%Create the equations
eq=BlochEqns(H,L);


%----------------------------TIME EVOLUTION-------------------------------%

%Before solving ODEs, display variables needed to be substituted
eq.necessaryVariables();

%Define initial conditions
IC=zeros(4);

IC(1,1)=1/10;
IC(2,2)=1/2;
IC(3,3)=2/5;

%Define numerical values for necessary substitutions. Chosen here to elucidate effects of dark and bright states
gamma=1;
Rabi_L=1;
Rabi_m=1;
delta_L=0;
delta_m=0;
Delta=0;

a_val=1/sqrt(2); 
b_val=1/sqrt(2);
A_val=1/sqrt(2);
B_val=-1/sqrt(2);

t_initial=0;
t_final=50;


%Solve the master equation
eq.evolve(t_initial,t_final,IC,[A_val,B_val,Delta,Rabi_L,Rabi_m,a_val,b_val,delta_L,delta_m,gamma]);

%And plot it
eq.plotEvolution();


%Let's check it in the bright and dark state basis
u2=subs(u2,[A,B],[A_val,B_val]);
u3=subs(u3,[A,B],[A_val,B_val]);
eq.changeBasis(u1,u2,u3,u4);

%And plot it in new basis. As you can see, the excited state population is
%0 after t~20, because the bright state is completely depopulated. The
%remaining poulation keeps on oscillating between the auxiliary state and
%the dark state.
figure
plot(eq.evTime,squeeze(eq.evolutionTr(1,1,:)))
hold on
plot(eq.evTime,squeeze(eq.evolutionTr(2,2,:)))
hold on
plot(eq.evTime,squeeze(eq.evolutionTr(3,3,:)))
hold on
plot(eq.evTime,squeeze(eq.evolutionTr(4,4,:)))
legend('Auxiliary state','Bright state','Dark state','Excited state');
drawnow








