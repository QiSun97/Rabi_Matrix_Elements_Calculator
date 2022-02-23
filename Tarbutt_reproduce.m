clear all;
import additional.*
import bloch.*
import BaH.* %Small module that includes state generation, dipole transition matrix elements calculation and Zeeman effect calculation in Barium Hydride
warning('off','all');

%***********************PHYSICAL CONSTANTS********************************%

%Frequency units used in this calculation are GHz, and time units of ns.

hbar=1.054*10^(-34); %[Js]
k_b=1.381*10^(-23); % [J/K]
c=299792000; %[m/s]
eps_0=8.854*10^(-12); %[F/m]
Gamma=0.964/58; %1/Lifetime of the excited state [GHz]
a0=5.29*10^(-11); % [m]
q_e=1.602*10^(-19); % [C]
r_expval=1.632*a0; % [m]  %need to ask
g_factor=-0.5/1000; %[GHz/G]  %need to ask
Bohr_mag=1.39962449/1000; %[GHz/G]

%*************************STATE GENERATION********************************%

[StG,StE,nground,nexcite]=generateSimpleStates(2,1); %For reference, check +BaH/generateStates.m
n=nground+nexcite; %Total number of states

%*************************TRANSITION DIPOLES******************************%

disp('Transition matrix')
rabi_matrix=zeros(n,n,3);
ind=0;
for p=[-1,0,1] %Three basic polarizations: -1 - right circular, +1 - left circular, 0 - pi
    ind=ind+1;
    for f=1:nexcite     %f - excited electronic states
        for i=1:nground %i - ground electronic states
            rabi_matrix(i,nground+f,ind)=dipoleTransitionMatrixElement(StG(i,:),StE(f,:),p); %For reference, +BaH/dipoleTransitionMatrixElement.m
            rabi_matrix(nground+f,i,ind)=dipoleTransitionMatrixElement(StE(f,:),StG(i,:),p);
        end
    end
end

%*************************ZEEMAN EFFECT***********************************%

disp('Zeeman effect matrix')
zeeman_matrix=sym(zeros(n,n));

syms g_g g_e B real; %Symbolic variables for g-factor and magnetic field
for i=1:nground     %add Zeeman shift to ground states
    zeeman_matrix(i,i)=g_g*B*StG(i,end);
end
for i=1:nexcite     %add Zeeman shift to excited states
    j=i+nground;
    zeeman_matrix(j,j)=g_e*B*StE(i,end);
end

%We calculate the Zeeman effect directly only for the ground state, which
%is in Hund's case (b) and for which we found appropriate formulas
% for i=1:nground     %f - final states
%     for f=1:nground %i - initial states
%         zeeman_matrix(i,f)=zeemanElement(StG(i,:),StG(f,:)); %For reference, +BaH/zeemanElement.m
%     end
% end

%*************************BRANCHING RATIOS********************************%

disp('Branching ratios')
BR=zeros(n);

%Transition strengths
transition_strengths=zeros(n);
for i=1:n
    for f=1:n
        for p=1:3
            transition_strengths(i,f)=transition_strengths(i,f)+rabi_matrix(i,f,p)^2;
        end
    end
end

%Sums of transition strengths for a given initial state 'i'
for i=1:n
    sums=0;
    for f=1:n
        sums=sums+transition_strengths(i,f);
    end
    for f=1:n
        BR(i,f)=transition_strengths(i,f)/sums; %(rotational) branching ratio
    end
end

%Initial states don't decay, so we remove those terms. Otherwise, they
%would simply indicate fractional transition strengths.
for i=1:nground
    BR(i,:)=0;
end

disp(BR(nground+1:end,:)) %Rows indicate branching ratios of decays from a single excited state (columns are ground states)

%****************************DISSIPATOR***********************************%

disp('Dissipator')
L=Dissipator(n);

syms G real;
assume(G,'positive')

%All excited states are assumed to have the same total decay rate G.
DR=zeros(1,nground);
for i=1:nexcite
    DR=[DR,G];
end
%We use branching ratios table to generate the whole dissipator, instead of
%adding decays one-by-one. That's also why we needed DR vector.
L.fromBranching(BR,DR);

%%
%****************************HAMILTONIAN**********************************%

%Symbolic variables
syms w_0 w_e real; %energies, which eventually will be replaced by detuning
syms w_L real; %light frequencies
syms W_L real; %Rabi frequencies.
syms d_L real; %detuning frequencies.
syms x v k t real; %additional variables

%Hamiltonian
disp('Hamiltonian')
H=Hamiltonian(n);

%States are added in following order:
%XSigma %1-5 - J=3/2, F=2, m=-2,-1,0,1,2
%BSigma %6-8 - J=1/2, F=1, m=-1,0,1
H.addEnergies([w_0,w_0,w_0,w_0,w_0,w_e,w_e,w_e]);

%add coupling by prefactor * E field modulation * rabi rate * rabi matrix element
phi = pi/4;
for i=1:nground
    for f=nground+1:n
        if rabi_matrix(i,f,1)~=0
            prefactor=(-1)^(-1)*sqrt(2);
            H.addPolyCoupling(i,f,prefactor*W_L*rabi_matrix(i,f,1)*cos(k*x-phi/2),w_L);
        end
        if rabi_matrix(i,f,3)~=0
            prefactor=-(-1)^(1)*sqrt(2);
            H.addPolyCoupling(i,f,prefactor*W_L*rabi_matrix(i,f,3)*cos(k*x+phi/2),w_L);
        end
        %         if rabi_matrix(i,f,2)~=0
        %             prefactor=(-1)^(0)*sqrt(2);
        %             H.addPolyCoupling(i,f,prefactor*W_L*rabi_matrix(i,f,2)*cos(k*x-phi/2),w_L);
        %         end
    end
end

H.hamiltonian=H.hamiltonian+zeeman_matrix;
H.defineEnergyDetuning(w_0,w_e,d_L,w_L);
H.defineZero(w_0);
H.unitaryTransformation();

for i=1:n
    for f=1:n
        if i~=f
            H.transformed(i,f)=simplify(expand(H.transformed(i,f)),'Steps',1000);
        end
    end
end

H.takeGradient(x);
H.transformed=simplify(subs(H.transformed,x,v*t));
DH=conj(simplify(subs(H.hamGradient,x,v*t)));

disp(H.transformed)
disp(DH)

%*************************INITIAL CONDITIONS******************************%

IC=zeros(n);
for i=1:nground
    IC(i,i)=1/nground;
end

%%
%**************************NUMERICAL VALUES*******************************%

%Decay rates
Gamma_BS=1/33; %1/Lifetime of the excited state [GHz]

%k-values
k_bs=2*pi/695; %1/nm

%laser freq detuning (which already replaced laser frequency and state energy)
det_laser = -2.5*Gamma_BS;

%Magnetic field
g_eff_G=0.56*Bohr_mag;
g_eff_E=0*Bohr_mag;
B_field=2*pi*12; %[G]

%Average Rabi rates [GHz]
% Rabi_L=-det_laser; %free parameter, should be calculated by I.
Rabi_L=sqrt(nexcite/2)*Gamma_BS; %free parameter, should be calculated by I.

t_equil=-2000/Gamma_BS; %[ns]

%*************************MASTER EQUATION*********************************%

disp('Optical Bloch Equations')
Eq=BlochEqns(H,L);
Eq.necessaryVariables();
disp(symvar(DH))
DHs=vpa(simplify(subs(DH,[W_L, k],[Rabi_L, k_bs])));

%%
%*************************VELOCITY PROFILE********************************%
% VelocitiesL=[-20:5:-10, -9:1:-1, -0.9:0.1:-0.1, -0.09:0.01:-0.01, 0];
% VelocitiesL=[-9:1:-1, -0.9:0.3:-0.1, -0.01, -0.005];
VelocitiesL=[0];
% VelocitiesL=[-0.01:0.001:0];
% VelocitiesL=[-20:5:-10, -9:1:-1, -0.9:0.1:-0.1, -0.095:0.005:-0.01, -0.0095:0.0005:0];

% VelocitiesR=-flip(VelocitiesL);
VelocitiesR=[];
Velocities=Gamma_BS/k_bs*[VelocitiesL,VelocitiesR];
nloop = length(Velocities);
%% Here comes the actual evaluation!
Eq_parfor=[];
for temp=1:nloop
    Eq_parfor=[Eq_parfor;copy(Eq)];
end

ppm = ParforProgressbar(nloop,'showWorkerProgress', true, 'progressBarUpdatePeriod', 1.5);

parfor vel_index=1:nloop
    tic
    velocity=Velocities(vel_index);
    t_end = 2000/Gamma_BS;
    Eq_parfor(vel_index).evolve(t_equil,t_end,IC,[B_field, Gamma_BS, Rabi_L, det_laser, g_eff_E, g_eff_G, k_bs, velocity]);
    Eq_parfor(vel_index).steadyState = [];
    toc
    ppm.increment();
end
delete(ppm);

%% plot excited state population vs time
% checklist = [1:2];
checklist = 1:nloop;
ncheck = length(checklist);

figure;
for ii=1:ncheck
    i = checklist(ii);
    eq = Eq_parfor(i);
    vel = Velocities(i);
    % cal & plot total excited state population
    ex_pop=real(squeeze(eq.evolution(6,6,:))...
        +squeeze(eq.evolution(7,7,:))+squeeze(eq.evolution(8,8,:)));
    plot(eq.evTime*Gamma_BS,ex_pop,'DisplayName',string(vel/(Gamma_BS/k_bs))+' * (Gamma/k)')
    hold on;
end
xlabel('time (/Gamma)')
ylabel('excited state total population')
legend
% xlim([0, eq.evTime(end)])

% return
%%
% checklist = [7];
checklist = 1:nloop;
ncheck = length(checklist);

figure;
tiledlayout('flow')

for ii=1:ncheck
    i = checklist(ii);
    eq = Eq_parfor(i);
    vel = Velocities(i);

    % perform FFT to total excited state population, to extract period
    [fft_init,step,force_endtime,fft_T]=FFTforce(nground, n, eq);
    fprintf('eigenfrequency is %.5f \Gamma \n',1/(fft_T * (eq.evTime(fft_init+1)-eq.evTime(fft_init)) * Gamma_BS))
    fprintf('period is %.2f, step is %.1f, start is %.1f, end is %.1f, length is %.1f \n',fft_T,step,fft_init,force_endtime,force_endtime-fft_init)
    
    Forces = [];
    tic
    for index=fft_init:step:force_endtime
        gradH=double(subs(DHs,[t,v],[eq.evTime(index),vel]));
        force=-real(trace(squeeze(eq.evolution(:,:,index))*gradH));
        Forces=[Forces,force];
    end
    toc
    tot_force=double(trapz(eq.evTime(fft_init:step:force_endtime),Forces)/(eq.evTime(force_endtime)-eq.evTime(fft_init)));
%     tot_force=double(mean(Forces));
    Forces_array(ii)=tot_force;
    fprintf('Velocity %.4f m/s\n',vel);
    fprintf('F=%.4f hbarkG*1E-3 \n',tot_force/(k_bs*Gamma_BS*1E-3));

    nexttile
    plot(eq.evTime(fft_init:step:force_endtime),Forces/(k_bs*Gamma_BS*1E-3));
    ylabel({'force at vel';num2str(vel)});
end

figure;
plot(Velocities(checklist)/(Gamma_BS/k_bs), Forces_array/(k_bs*Gamma_BS*1E-3));
xlabel('velocity (\Gamma/k)');
ylabel('force (k*\Gamma*1E-3)')
% xlim([-20, 20])
set(gcf, 'position', [100, 100, 300, 200])

return
%%
% checklist = [7];
checklist = 1:nloop;
ncheck = length(checklist);

for ii=1:ncheck
    i = checklist(ii);
    eq = Eq_parfor(i);
    vel = Velocities(i);
    
    % perform FFT to total excited state population, to extract period
    ex_pop=real(squeeze(eq.evolution(6,6,:))...
        +squeeze(eq.evolution(7,7,:))+squeeze(eq.evolution(8,8,:)));
    fft_init=find(eq.evTime(:)>0,1);
    fft_endtime = length(eq.evTime(:));
    fft_endtime = fft_endtime-1+mod(fft_endtime-fft_init,2);
    fft_data = ex_pop(fft_init:fft_endtime);
    fft_datalength = fft_endtime-fft_init+1;
    fft_Y = abs(fft(fft_data));
    fft_cutoff = 2;
    fft_Y = fft_Y(fft_cutoff+1:(fft_datalength/2)+1);
    fft_freq = (fft_cutoff:(fft_datalength)/2)/fft_datalength;
    [temp,peak_freq] = max(fft_Y);
    fft_T = 1/(fft_freq(peak_freq));
    
    if fft_T>2000
        step = floor(fft_T/2000)+1;
        force_endtime = floor(fft_T)-1+fft_init;
    elseif fft_T<600
        step = 1;
        force_endtime = floor(floor(600/fft_T)*fft_T)-1+fft_init;
    else
        step = 1;
        force_endtime = floor(fft_T)-1+fft_init;
    end
    if force_endtime>=fft_endtime
        disp('something is wrong, your force end time is too long')
        force_endtime = fft_endtime;
    end
    fprintf('eigenfrequency is %.5f \Gamma \n',1/(fft_T * (eq.evTime(2)-eq.evTime(1)) * Gamma_BS))
    fprintf('period is %.2f, step is %.1f, start is %.1f, end is %.1f, length is %.1f \n',fft_T,step,fft_init,force_endtime,force_endtime-fft_init)
    Forces = [];
    tic
    for index=fft_init:step:force_endtime
        gradH=double(subs(DHs,[t,v],[eq.evTime(index),vel]));
        force=-real(trace(squeeze(eq.evolution(:,:,index))*gradH));
        Forces=[Forces,force];
    end
    tot_force=double(trapz(eq.evTime(fft_init:step:index),Forces)/(eq.evTime(index)-eq.evTime(fft_init)));
%     tot_force=double(mean(Forces));
    Forces_array(ii)=tot_force;
    toc
    fprintf('Velocity %.4f m/s\n',vel);
    fprintf('F=%.4f hbarkG*1E-3 \n',tot_force/(k_bs*Gamma_BS*1E-3));
end

figure;
plot(Velocities(checklist)/(Gamma_BS/k_bs), Forces_array/(k_bs*Gamma_BS*1E-3));
xlabel('velocity (\Gamma/k)');
ylabel('force (k*\Gamma*1E-3)')
% xlim([-20, 20])
set(gcf, 'position', [100, 100, 300, 200])

% return
%% calculate force vs vel and time
% checklist = 1:10;
checklist = 1:nloop;
ncheck = length(checklist);

% DHs_parfor=[];
% for temp=1:ncheck
%     DHs_parfor=[DHs_parfor;DHs];
% end

Forces_array = zeros(ncheck,1);
% ppm = ParforProgressbar(ncheck,'showWorkerProgress', true, 'progressBarUpdatePeriod', 1.5);

for ii=1:ncheck
    tic
    i = checklist(ii);
    eq = Eq_parfor(i);
    vel = Velocities(i);
    
    i0=find(eq.evTime(:)>0,1);
    step = 1;
    %     lengthtime = length(eq.evTime(:));
    lengthtime = i0+200;
    
    Forces = [];
    for index=i0:step:lengthtime
        gradH=double(subs(DHs,[t,v],[eq.evTime(index),vel]));
        force=-real(trace(squeeze(eq.evolution(:,:,index))*gradH));
        Forces=[Forces,force];
    end
    %     tot_force=double(trapz(eq.evTime(i0:step:lengthtime),Forces)/(eq.evTime(index)-eq.evTime(i0)));
    tot_force=double(mean(Forces));
    Forces_array(ii)=tot_force;
    toc
    %     ppm.increment(s);
    fprintf('Velocity %.4f m/s\n',vel);
    fprintf('F=%.4f hbarkG*1E-3 \n',tot_force/(k_bs*Gamma_BS*1E-3));
end

figure;
plot(Velocities(checklist)/(Gamma_BS/k_bs), Forces_array/(k_bs*Gamma_BS*1E-3));
xlabel('velocity / (Gamma/k)');
ylabel('force / (k*Gamma*1E-3)')
xlim([-20, 20])
set(gcf, 'position', [100, 100, 300, 200])

return
%% save arrays
save('0617data3.mat','Velocities','Gamma_BS','k_bs',...
    'Forces_array','phi','H','Eq');

%% plot using saved data
checklist=1:length(Velocities);
% figure;
plot([Velocities(checklist)/(Gamma_BS/k_bs),...
    flip(-Velocities(checklist))/(Gamma_BS/k_bs)],...
    [-Forces_array/(k_bs*Gamma_BS*1E-3),flip(Forces_array)/(k_bs*Gamma_BS*1E-3)],...
    'DisplayName','phi = 0');
xlabel('velocity / (Gamma/k)');
ylabel('force / (k*Gamma*1E-3)');
legend;
% xlim([-5, +5]);
% ylim([-7, 7]);
set(gcf, 'position', [100, 100, 500, 400])
hold on

%% some test code
% checklist = [8];
checklist = 1:nloop;
ncheck = length(checklist);

for ii=1:ncheck
    i = checklist(ii);
    eq = Eq_parfor(i);
    vel = Velocities(i);
    disp(['vel is ',num2str(vel)])
    
    % perform FFT to total excited state population, to extract period
    ex_pop=real(squeeze(eq.evolution(6,6,:))...
        +squeeze(eq.evolution(7,7,:))+squeeze(eq.evolution(8,8,:)));
    fft_init=find(eq.evTime(:)>0,1);
    fft_endtime = length(eq.evTime(:));
    fft_endtime = fft_endtime-1+mod(fft_endtime-fft_init,2);
    fft_data = ex_pop(fft_init:fft_endtime);
    fft_datalength = fft_endtime-fft_init+1;
    fft_Y = abs(fft(fft_data));
    fft_cutoff = 2;
    fft_Y = fft_Y(fft_cutoff+1:(fft_datalength/2)+1);
    fft_freq = (fft_cutoff:(fft_datalength)/2)/fft_datalength;
    [temp,peak_freq] = max(fft_Y);
    fft_T = 1/(fft_freq(peak_freq));
    
    if fft_T>600
        step = int(fft_T/800)+1;
        force_endtime = floor(fft_T)-1+fft_init;
    elseif fft_T<200
        step = 1;
        force_endtime = floor(floor(200/fft_T)*fft_T)-1+fft_init;
    else
        step = 1;
        force_endtime = floor(fft_T)-1+fft_init;
    end
    
    fprintf('period is %.2f, step is %.1f, start is %.1f, end is %.1f, length is %.1f',fft_T,step,fft_init,force_endtime,force_endtime-fft_init)
    Forces = [];
    for index=fft_init:step:force_endtime
        gradH=double(subs(DHs,[t,v],[eq.evTime(index),vel]));
        force=-real(trace(squeeze(eq.evolution(:,:,index))*gradH));
        Forces=[Forces,force];
    end
    tot_force=double(mean(Forces));
    Forces_array(ii)=tot_force;
    toc
    fprintf('Velocity %.4f m/s\n',vel);
    fprintf('F=%.4f hbarkG*1E-3 \n',tot_force/(k_bs*Gamma_BS*1E-3));
end

figure;
plot(Velocities(checklist)/(Gamma_BS/k_bs), Forces_array/(k_bs*Gamma_BS*1E-3));
xlabel('velocity / (Gamma/k)');
ylabel('force / (k*Gamma*1E-3)')
xlim([-20, 20])
set(gcf, 'position', [100, 100, 300, 200])

%% some test code
checklist = [4, 29, 30];
% checklist = 1:nloop;
ncheck = length(checklist);

for ii=1:ncheck
    i = checklist(ii);
    eq = Eq_parfor(i);
    vel = Velocities(i);
    disp(vel)
    
    % cal & plot total excited state population
    ex_pop=real(squeeze(eq.evolution(6,6,:))...
        +squeeze(eq.evolution(7,7,:))+squeeze(eq.evolution(8,8,:)));
    %     figure;
    %     plot(eq.evTime*Gamma_BS,ex_pop,'DisplayName',string(vel/(Gamma_BS/k_bs))+' * (Gamma/k)')
    fft_init=find(eq.evTime(:)>0,1);
    fft_endtime = length(eq.evTime(:));
    %     fft_endtime = fft_endtime-1+mod(fft_endtime-fft_init,2);
    fft_endtime = fft_init+599;
    fft_data = ex_pop(fft_init:fft_endtime);
    fft_datalength = fft_endtime-fft_init+1;
    fft_Y = abs(fft(fft_data));
    fft_Y = fft_Y(1:(fft_datalength/2)+1);
    fft_freq = (0:(fft_datalength)/2)/fft_datalength;
    
    figure;
    subplot(2,2,1)
    plot(fft_freq(1:end),fft_Y(1:end));
    title('pop FFT')
    
    Forces = [];
    for index=fft_init:fft_endtime
        gradH=double(subs(DHs,[t,v],[eq.evTime(index),vel]));
        force=-real(trace(squeeze(eq.evolution(:,:,index))*gradH));
        Forces=[Forces,force];
    end
    
    forcefft_freq = (0:(fft_datalength)/2)/fft_datalength;
    forcefft_Y = abs(fft(Forces));
    forcefft_Y = forcefft_Y(1:fft_datalength/2+1);
    
    subplot(2,2,2)
    plot(forcefft_freq,forcefft_Y);
    title('force FFT')
    
    subplot(2,2,3)
    plot(eq.evTime(fft_init:fft_endtime),Forces)
    title('pop time')
    xlabel('time (/Gamma)')
    
    subplot(2,2,4)
    plot(eq.evTime(fft_init:fft_endtime),ex_pop(fft_init:fft_endtime))
    title('force time')
    xlabel('time (/Gamma)')
    %     hold on;
end
% ylabel('excited state total population')
% legend

%% plot excited state population vs time
checklist = [30];
% checklist = 1:nloop;
ncheck = length(checklist);
figure;
for ii=1:ncheck
    i = checklist(ii);
    eq = Eq_parfor(i);
    vel = Velocities(i);
    % cal & plot total excited state population
    ex_pop=real(squeeze(eq.evolution(6,6,:))...
        +squeeze(eq.evolution(7,7,:))+squeeze(eq.evolution(8,8,:)));
    plot(eq.evTime*Gamma_BS,ex_pop,'DisplayName',string(vel/(Gamma_BS/k_bs))+' * (Gamma/k)')
    hold on;
end
xlabel('time (/Gamma)')
ylabel('excited state total population')
legend
%%
% i0=find(Eq_inpar.evTime(:)>0,1);
%
% av_e=Eq_inpar.evolution(nx+1,nx+1,i0:end);
% for ii=nx+2:n
%     av_e=av_e+Eq.evolution(ii,ii,i0:end);
% end
% av_e=mean(av_e);
%
% for i=i0:step:length(Eq.evTime(:))
%     gradH=double(subs(DHs,[t,v],[Eq.evTime(i),velocity]));
%     force=-real(trace(squeeze(Eq.evolution(:,:,i))*gradH));
%     Forces=[Forces,force];
% end
%
% tot_force=double(trapz(Eq.evTime(i0:step:i),Forces)/(Eq.evTime(i)-Eq.evTime(i0)));
% tot_force=tot_force/(k_api*Gamma_APi*0.5);
% Forces_array_extended_more(vel_index)=tot_force;
% fprintf('Velocity %.2f m/s\n',velocity);
% fprintf('p_ee=%.4f \n',av_e);
% fprintf('F=%.4f hkG/2 \n',tot_force);
% fprintf('Done %.2f%% \n',u/(31)*100);
%
% Eq.plotEvolution();
%
%
% Fs=[Fs,tot_force];
%
% clf
% plot(-10:2:velocity,Fs)
% xlim([-10 10])
% ylim([-40 40])
% drawnow