clc
close all
clear all

%% geometria

D_p=1e-2;%[m],  diametro plasma
D_nozz=5e-2; %[m] , diametro ugelli
A_nozz=D_nozz^2*pi/4; % [m^2], superficie ugelli



% camera centrale
R_c=0.13/2; %[m], raggio
L_c=0.16; %[m], lunghezza
As_c=2*pi*R_c*L_c; %[m^2], superficie laterale
s_c=0.01; %[m] , spessore cps (ipotizato)
V_c_wall = pi*(R_c^2-(R_c-s_c)^2)*L_c; % [m^3], volume cps
V_c_box=pi*(R_c-s_c)^2*L_c; %[m^3], volume interno box
A_inner_c=D_p*L_c*pi; % [m^2], area di interazione del plasma

% camera upstream (sinistra)
R_sx=0.13/2; %[m], raggio
L_sx=0.08; %[m], lunghezza
As_sx=2*pi*R_sx*L_sx; %[m^2], superficie laterale
s_sx=0.01; %[m] , spessore cps (ipotizzato)
V_sx_wall = pi*(R_sx^2-(R_sx-s_sx)^2)*L_sx; % [m^3], volume cps
V_sx_box=pi*(R_sx-s_sx)^2*L_sx; %[m^3], volume interno box
A_inner_s=pi*D_p*L_sx; % [m^2], area di interazione del plasma

% camera downstream (destra)
R_dx=0.105/2; %[m], raggio
L_dx=0.16; %[m], lunghezza
As_dx=2*pi*R_dx*L_dx; %[m^2], superficie laterale
s_dx=0.01; %[m] spessore cps (ipotizzato)
V_dx_wall = pi*(R_dx^2-(R_dx-s_dx)^2)*L_dx; % [m^3], volume cps
V_dx_box=pi*R_dx^2*L_dx; %[m^3], volume interno box
A_inner_d=pi*D_p*L_dx; %[m^2], area di interazione del plasma

% target
D_target= 0.05; %[m]
s_cu= 4e-3; %[m]
s_tung= 4e-3; %[m]
s_target=s_cu+s_tung; %[m]
As_target= pi*D_target^2/4; %[m^2]

%% cps

% Proprietà Litio
% T in K
rho_li = @(T) (278.5-0.04657.*T+274.6.*(1-T/3500)^0.467); % [kg/m^3] 
cp_li = @(T) (4754-0.925.*T+2.91*10^-4.*T.^2); %[J/(kg*K)]
k_li= @(T) (22.28+0.0500.*T-1.243*1e-5.*T^2); %[W/(m*K)]

% Proprietà Tungsteno
% T in K
rho_tung=@(T) (1e3.*(19.3027-2.3786e-4.*(T-273.15)-2.2448e-8.*(T-273.15).^2)); % [kg/m^3] 
cp_tung=@(T) (128.308+3.2797e-2.*(T-273.15)-3.4097e-6.*(T-273.15).^2); % [J/(kg*K)]
k_tung=@(T) (174.9274-0.1067.*(T-273.15)+5.0067e-5.*(T-273.15).^2-7.8349e-9.*(T-273.15).^3); % [W/(m*K)] 

%CPS
f=0.5;
k_CPS=@(T) f*k_li(T)+(1-f)*k_tung(T); %[W/m/K]
rho_CPS=@(T) f*rho_li(T)+(1-f)*rho_tung(T);  %[kg/m^3]
cp_CPS=@(T) f*cp_li(T)+(1-f)*cp_tung(T); %[J/kg/K]

% Proprietà Rame
MMcu=63.54e-3; %[Kg/mol]
rho_cu= @(T) (1e3.*(8.81-4.28e-4.*(T-300)-6.12e-8.*(T-300).^2-2.77e-11.*(T-300).^3)); % [kg/m^3] 
cp_cu= @(T) ((24.27+1.23e-2.*(T-300)-2.05e-5.*(T-300).^2+1.53e-8.*(T-300).^3-2.88e-12.*(T-300).^4)./MMcu); %[J/(kg*K)]
k_cu= @(T) (11.627-4.2781e-2.*(T-300)+9.2403e-5.*(T-300).^2-9.6577e-8.*(T-300).^3+3.686e-11.*(T-300).^4); %[W/(m*K)]

%costanti
K=physconst('Boltzmann'); %J/K costante Boltzmann
m=1.1525801e-26; %Kg (massa di un atomo di litio in Kg)
M=6.941e-3; %Kg/mol (massa molare di litio)
f_redep=0.9; % coefficiente tasso di evaporazione
eta=0.75; 
Na=6.022*1e23; %[mol^-1] Costante di Avogadro

%% funzioni
enthalpy_li_int= @(T) (5*K.*T/2); %[J/#], entalpia litio per atomo
p_v= @(T) (exp(26.89-18880./T-0.4942.*log(T))); %[Pa] T in K, tensione di vapore 
ev = @(T,A) (A*f_redep*eta*p_v(T)/sqrt(2*pi*m*K*T)); % [#/s], funzione di evaporazione
cond= @(Tvap,n_i,A) (A.*n_i.*sqrt(8.*K.*Tvap./(pi*m))./4); %[#/s] , funzione di condensazione

%% Heater
q_h=2000;% [W] 

%% Irraggiamento con esterno
sigma=5.67*1e-8; %costante di Stefan-Boltzmann [W/(m^2*K^4)]
epsilon=1; % emissività del corpo nero (ipotesi)  
Tp=20+273.15; %[K] %temperatura parete del laboratorio

%% Acqua target
T_water=293.15; %[K]
pressure=5; %[MPa]
T_sat=536.95; %[K]
rho_water=1000; %[Kg/m^3]
v_water=12; %[m/s]
Pr_water=1.20;
k_water=688e-3; %[W/(m*K)]
D_water=0.01; %[m]
mu_water=8.9e-4; %[Pa*s]
mu_d_water=mu_water/rho_water;
Re_water=rho_water*v_water*D_water/mu_d_water;
Nu_water=0.023*Re_water^(4/5)*Pr_water^0.4;
hh_w=Nu_water*k_water/D_water;

%% Condizioni iniziali
% Temperature pareti delle camere
Twall_c=293.15; %[K]
Twall_dx=293.15; %[K]
Twall_sx=293.15; %[K]

% Flusso volumico
q_F_c=0/V_c_wall; % [W/m^3]
q_F_dx=0/V_dx_wall; % [W/m^3]
q_F_sx=0/V_sx_wall; % [W/m^3]

% Concentrazione di litio
Ncold=0; % [#/m^3]
Nsold=0; % [#/m^3]
Ndold=0; % [#/m^3]

%Litio in CPS
MMc=100e-3*Na/M;  % [kg]

% Temperatura del vapore
Tvcold=293.15; % [K]
Tvdold=293.15; % [K]
Tvsold=293.15; % [K]

% Concentrazione di litio iniziale nelle camere [sx c dx]
n_pre_V = [0 0 0]; 

% [Temperatura elettronica target; Temperatura el. upstream; Flusso areico
% al target] iniziale
sol_pre_v = [20; 24; 4e7]; % [eV; eV; W/m^2]

% parametro per la soluzione numerica del sistema di Lengyel
steps = 1; 

% iterazioni massime
i_max = 1000;

% preallocazione del vettore tempo
t = zeros(i_max,1);

%% Discretizzazione
dt=5; % [s]
Tempi=[3000 3700 4200 4660];
jj=1;
i=1;

%discretizzazione target
dx=1e-5;
x_tung=(0:dx:s_tung)';
x_cu=(s_tung+dx:dx:s_target)';
x_target=[x_tung;x_cu];
N_tung=length(x_tung);
N_cu=length(x_cu);
N_target=length(x_target);

% Temperatura target iniziale
Told_tung=293.15*ones(N_tung,1);
Told_cu=293.15*ones(N_cu,1);
Told_target=[Told_tung;Told_cu];

%per i plot
Tc(i)=Tvcold;
Ts(i)=Tvsold;
Td(i)=Tvdold;
MM(i)=MMc*M/Na;
lost_li(i)=0;
lost_li_complessivo = zeros(i_max,1);

%% Ciclo iterativo

% dizionario concentrazioni di litio nelle 3 camere
dizionario_nz = [1e20 1e20 1e20; 1e22 1e22 1e22; 2e22 2e22 2e22; 5e22 5e22 5e22; 1e23 1e23 1e23; 1e25 1e25 1e25]; 
% dizionario soluzione del sistema per le concentrazioni di litio del
% precedente dizionario
dizionario_sol = [1.38,10.47,4.8e+06;0.21 14 2.5e6; 1.39e-4  13  6.5e4; 7.8e-6  11  1.4e4 ; 6.4e-6 10 7.7e3;1.1e-4 2.3 1.3e3];

% Inizializzazione flags
failedFlag = 0;
dictionaryFlag = 0;

% Prealloco per conservare le soluzioni nelle varie iterazioni 
prima_colonna_AA = zeros(i_max,3);
seconda_colonna_AA = zeros(i_max,3);
terza_colonna_AA = zeros(i_max,3);
registro_Fail = [];
registro_dict = [];

while  i<i_max
    i=i+1;
    t(i)=t(i-1)+dt;


    %% Parete
    g0=dt./(rho_CPS(Twall_sx(i-1)).*cp_CPS(Twall_sx(i-1)));
    g1=dt./(rho_CPS(Twall_sx(i-1)).*cp_CPS(Twall_sx(i-1)).*V_sx_wall);
    g2=dt*As_sx*sigma*epsilon./(rho_CPS(Twall_sx(i-1)).*cp_CPS(Twall_sx(i-1)).*V_sx_wall);
    
    j0=dt./(rho_CPS(Twall_dx(i-1)).*cp_CPS(Twall_dx(i-1)));
    j1=dt./(rho_CPS(Twall_dx(i-1)).*cp_CPS(Twall_dx(i-1)).*V_dx_wall);
    j2=dt*As_dx*sigma*epsilon./(rho_CPS(Twall_dx(i-1)).*cp_CPS(Twall_dx(i-1)).*V_dx_wall);
    
    b0=dt./(rho_CPS(Twall_c(i-1)).*cp_CPS(Twall_c(i-1))); 
    b1=dt./(rho_CPS(Twall_c(i-1)).*cp_CPS(Twall_c(i-1)).*V_c_wall); 
    b2=dt./(rho_CPS(Twall_c(i-1)).*cp_CPS(Twall_c(i-1)).*V_c_wall); 
    
    
     % controllo se è presente litio nella cps per l'evaporazione
    if MMc>0 
        if MMc-ev(Twall_c(i-1),As_c)*dt>0 
            
            ev_c=ev(Twall_c(i-1),As_c); % [#/s]
            
        else
            ev_c=MMc/dt;% [#/s]
           
        end
    else
        ev_c=0;
    end
    

    % Temperature delle pareti all'iterazione
    
    Twall_c(i)=b1*q_h+Twall_c(i-1)-b2*enthalpy_li_int(Twall_c(i-1))*ev_c + q_F_c*b0+b2*enthalpy_li_int(Tvcold)*cond(Tvcold,Ncold,As_c);   % [K]

    Twall_dx(i)=j2*Tp^4-j2.*Twall_dx(i-1)^4+Twall_dx(i-1)+ q_F_dx*j0+j1*enthalpy_li_int(Tvdold)*cond(Tvdold,Ndold,As_dx); % [K]

    Twall_sx(i)=g2*Tp^4-g2.*Twall_sx(i-1)^4+Twall_sx(i-1)+ q_F_sx*g0+g1*enthalpy_li_int(Tvsold)*cond(Tvsold,Nsold,As_sx); % [K]

    %% Vapore
    % controllo se è presente litio nella cps per l'evaporazione
    if MMc>0 
        if MMc-ev(Twall_c(i),As_c)*dt>0 
            evc=ev(Twall_c(i),As_c); % [#/s]
            MMc=MMc-evc*dt;% [#/s]
        else
            evc=MMc/dt;% [#/s]
            MMc=0;
        end
    else
        evc=0;
    end

    %camera centrale
    y=2*(evc*5*K*Twall_c(i)/2*dt/V_c_box+Ncold*5*K*Tvcold/2)/(5*K*(evc*dt/V_c_box+Ncold)); %Tc (K)
    x=(evc*dt/V_c_box+Ncold)/(1+dt*sqrt(y)*(As_c*sqrt(8*K/(pi*m))/4+A_inner_c*sqrt(8*K/(pi*m))/4+2*0.6288*A_nozz*sqrt(K/m))/V_c_box); %Nc [#/m^3]
    n_c = x;
    %camera sinistra
    j=2*(0.6288*A_nozz*x*sqrt(K*y/m)*5*K*y/2*dt/V_sx_box+Nsold*5*K*Tvsold/2)/(5*K*(Nsold+0.6288*A_nozz*x*sqrt(K*y/m)*dt/V_sx_box)); %Ts [K]
    z=(Nsold+0.6288*A_nozz*x*sqrt(K*y/m)*dt/V_sx_box)/(1+sqrt(j)*dt*(As_sx*sqrt(8*K/(pi*m))/4+0.6288*A_nozz*sqrt(K/m)+A_inner_s*sqrt(8*K/(pi*m))/4)/V_sx_box); %Ns [#/m^3]
    n_sx = z;
    %camera destra
    J=2*(0.6288*A_nozz*x*sqrt(K*y/m)*5*K*y/2*dt/V_dx_box+Ndold*5*K*Tvdold/2)/(5*K*(Ndold+0.6288*A_nozz*x*sqrt(K*y/m)*dt/V_dx_box)); %Td [K]
    Z=(Ndold+0.6288*A_nozz*x*sqrt(K*y/m)*dt/V_dx_box)/(1+sqrt(J)*dt*(As_dx*sqrt(8*K/(pi*m))/4+0.6288*A_nozz*sqrt(K/m)+A_inner_d*sqrt(8*K/(pi*m))/4)/V_dx_box); %Nd [#/m^3]
    n_dx = Z;
    Nc(i)=x;
    Ns(i)=z;
    Nd(i)=Z;
    Tc(i)=y;
    Ts(i)=j;
    Td(i)=J;
    Ncold=x;
    Nsold=z;
    Ndold=Z;
    Tvcold=y;
    Tvsold=j;
    Tvdold=J;

    % litio rimasto in cps e perso 
    MMc=MMc+As_c*x*sqrt(8*K*y/(pi*m))/4*dt+As_dx*Z*sqrt(8*K*J/(pi*m))/4*dt+As_sx*z*sqrt(8*K*j/(pi*m))/4*dt; % particelle
    lost_li(i)=0.6288*A_nozz*z*sqrt(K*j/m)*M/6.022e23+0.6288*A_nozz*Z*sqrt(K*J/m)*M/6.022e23+sqrt(8*K/(pi*m))*M/Na/4*(A_inner_c*x*sqrt(y)+A_inner_s*z*sqrt(j)+A_inner_d*Z*sqrt(J)); %[Kg/s]
    MM(i)=MMc*M/Na; %[Kg]
    
    
    %% Plasma

    % per il debugging
    % if i == 2500
    %     steps = 5;
    % end
    %     if sol_pre_v(1)<1
    %         oo = 'stop'
    %     end
    
    % risoluzione sistema con griglia a 200 nodi
    AA = plasma_fun(n_sx,n_c,n_dx,n_pre_V,sol_pre_v,steps,200);
    if AA(1,3) == 1
        % se non trova una soluzione prova con una griglia più raffinata
        fprintf("Fail 200\n")
        AA = plasma_fun(n_sx,n_c,n_dx,n_pre_V,sol_pre_v,steps,300);
        if AA(1,3) == 1
            fprintf("Fail 300\n")
            AA = plasma_fun(n_sx,n_c,n_dx,n_pre_V,sol_pre_v,steps,400);
            
            if AA(1,3) == 1
                fprintf("Fail 400\n")
                AA = plasma_fun(n_sx,n_c,n_dx,n_pre_V,sol_pre_v,steps,500);
                if AA(1,3) == 1
                    % se non trova una soluzione con 500 nodi, prova a
                    % prendere come guess il valore più vicino presente
                    % nel dizionario, il numero di nodi della griglia viene
                    % gestito nello stesso modo
                    fprintf("Fail 500\n")
                    distanze_soluzioni_mat = mean([n_sx n_c n_dx])-dizionario_nz;
                    distanze_soluzioni = mean(distanze_soluzioni_mat');
                    [dist_min,indexmin] = min(abs(distanze_soluzioni));
                    n_guess = dizionario_nz(indexmin,:);
                    sol_guess = dizionario_sol(indexmin,:);
    
                    AA = plasma_fun(n_sx,n_c,n_dx,n_guess,sol_guess',steps,200);
                    dictionaryFlag = dictionaryFlag + 1;
                    registro_dict = [registro_dict;i];
                    if AA(1,3) == 1
                        
                        fprintf("Fail 200 dict\n")
                        AA = plasma_fun(n_sx,n_c,n_dx,n_guess,sol_guess',steps,300);
                        if AA(1,3) == 1
                            fprintf("Fail 300 dict\n")
                            AA = plasma_fun(n_sx,n_c,n_dx,n_guess,sol_guess',steps,400);
                            if AA(1,3) == 1
                                fprintf("Fail 400 dict\n")
                                AA = plasma_fun(n_sx,n_c,n_dx,n_guess,sol_guess',steps,500);
                                if AA(1,3) == 1
                                    % se non trova una soluzione, rimane la
                                    % soluzione del nodo precedente,
                                    % aggiorno il registro e la flag
                                    fprintf("Fail 500\n")
                                    fprintf("Failed\n")
                                    failedFlag = failedFlag +1;
                                    registro_Fail = [registro_Fail; i];
                                    
                                end
                            end
                        end
                    end
                end
            end
         end
    end

% per il debugging
% if i == 931
%     oo = 'stop'
% end

    % Flussi volumici
    q_F_sx = (AA(1,1)*(AA(1,1)>0))/V_sx_wall;
    q_F_c = AA(2,1)*(AA(2,1)>0)/V_c_wall;
    q_F_dx = AA(3,1)*(AA(3,1)>0)/V_dx_wall;

    % aggiorno la soluzione del nodo per la risoluzione del sistema
    n_pre_V = [n_sx n_c n_dx];
    sol_pre_v = [AA(1,2)*(AA(1,2)>0);AA(2,2)*(AA(2,2)>0);AA(3,2)*(AA(3,2)>0)];
    
    % conservo i valori in dei registri
    prima_colonna_AA(i,:) = AA(:,1)';
    seconda_colonna_AA(i,:) = AA(:,2)';
    terza_colonna_AA(i,:) = AA(:,3)';

    i % stampo l'iterazione corrente

 %% Target
    %proprietà
    rho_t=rho_tung(Told_tung);
    rho_c=rho_cu(Told_cu);
    rho_target=[rho_t;rho_c];
    k_t=k_tung(Told_tung);
    k_c=k_cu(Told_cu);
    k_target=[k_t;k_c];
    cp_t=cp_tung(Told_tung);
    cp_c=cp_cu(Told_cu);
    cp_target=[cp_t;cp_c];
    
    a_t=k_target.*dt./(rho_target.*cp_target.*dx^2);
    
    % matrice AA
    sup=-a_t;
    main=1+2*a_t;
    sub=-a_t;
    
    aus=[[sub(2:end);0] main [0;sup(1:end-1)]];

    AA_target=spdiags(aus,-1:1,N_target,N_target);
    
    %%Condizioni al contorno
    AA_target(1,:)=0;
    AA_target(1,1)=-1;
    AA_target(1,2)=1;
    AA_target(end,:)=0;
    AA_target(end,end)=1+hh_w*dx/k_target(end);
    AA_target(end,end-1)=-1;

    AA_target(N_tung,N_tung+1)=-k_target(N_tung+1); 
    AA_target(N_tung,N_tung)=k_target(N_tung+1)+k_target(N_tung-1);
    AA_target(N_tung,N_tung-1)=-k_target(N_tung-1); 

    Tnoto=Told_target;
    qq_target=seconda_colonna_AA(i,3)*pi*D_p^2/4/(pi*D_target^2/4);
    Tnoto(1)=-dx/k_target(1)*qq_target;
    Tnoto(end)=hh_w*dx/k_target(end)*T_water;
    Tnoto(N_tung)=0;
    % soluzione
    Tnew_target=AA_target\Tnoto;
    Told_tung=Tnew_target(1:N_tung);
    Told_cu=Tnew_target(N_tung+1:N_target);
    Told_target=[Told_tung;Told_cu];
    
    if jj<=length(Tempi) && (i-1)*dt==Tempi(jj)
        figure(2)
        subplot(2,2,2)
        plot(x_target,Tnew_target)
        hold on
        grid on
        jj=jj+1;
    end

end

% salvo l'indice per cui la potenza al target è minima
[P_t_min,min_ind] = min(seconda_colonna_AA(2:end,3)); 


% salvo il vettore delle concentrazioni di litio e della soluzione del sistema all'istante trovato
min_nz = [Ns(min_ind),Nc(min_ind),Nd(min_ind)]; 
min_sol = seconda_colonna_AA(min_ind,:);
min_potenze_vettore = plasma__potenza_fun(min_nz(1),min_nz(2),min_nz(3),min_sol(1),min_sol(2),min_sol(3),200);

figure(1)
subplot(2,2,1)
plot(t,seconda_colonna_AA(:,3)*pi*D_p^2/4)
xlabel('t[s]')
ylabel('Potenza[W]')
title('Potenza al target')
grid minor

subplot(2,2,2)
plot(t,Nc)
hold on 
plot(t,Ns);
hold on 
plot(t,Nd)
ylabel('particelle/m^3');
xlabel('t [s]');
title('Densità litio ')
legend('Nc','Ns','Nd','Location','best')
grid minor

subplot(2,2,3)
plot(t,Twall_c)
hold on
plot(t,Twall_sx)
hold on
plot(t,Twall_dx)
ylabel('T [K]');
xlabel('t [s]');
title('Temperatura delle pareti ')
legend('Twall c','Twall sx','Twall dx','Location','best')
grid minor

subplot(2,2,4)
plot(t,Tc)
hold on
plot(t,Ts)
hold on
plot(t,Td)
ylabel('T [K]');
xlabel('t [s]');
title('Temperatura del vapore ')
legend('Tvap c','Tvap sx','Tvap dx','Location','best')
grid minor

figure(2)
subplot(2,2,1)
plot(linspace(0.52,0,length(min_potenze_vettore)),min_potenze_vettore*pi*D_p^2/4)
grid minor
xlabel('z (m)')
ylabel('Potenza (W)')
title('Distribuzione della potenza del plasma')
tempo_min = min_ind*dt;
subtitle_var = {'all''istante nel quale il plasma viene frenato maggiormente,',strcat('ovvero per t = ',num2str(tempo_min), ' s')};
subtitle(subtitle_var)

subplot(2,2,2)
xlabel('x[m]')
ylabel('T[K]')
title('Temperatura del target nei vari istanti di tempo')
legend('3000s', '3700s', '4200s','4660s')
grid minor

subplot(2,2,3)
plot(t,MM)
xlabel('t [s]')
ylabel('[Kg]')
title('Litio in CPS ')
grid minor

subplot(2,2,4)
plot(t,lost_li)
ylabel('[Kg/s]')
xlabel('t[s]')
title('portata litio perso')
grid minor


