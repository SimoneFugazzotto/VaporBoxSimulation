function AA = plasma_fun(nz1,nz2,nz3,nzPreV,solPreV,steps,nodi)
% accetta in m^-3, sputa in W
% prima colonna: potenza irraggiata nella camera upstream; pot irr camera
% centrale; pot irr camera targ [W]
% seconda colonna: temperatura target [eV]; temperatura upstream [eV]; flusso target [W/m^2]
% terza colonna: iter>maxIter; iter; maxIter

%% dati di input
flagShow = 0;
ne_u = 1e20; % m^-3
tau = 1e-3; % s
%qqu = 5e6; % W/m^2
ke_par = 2390; % W/m/eV^7/2
%ke_par = 2000;
AreaPlasma = pi*0.01^2/4;
%Pot = qqu*AreaPlasma;
Pot = 3300; % W
qqu = Pot/AreaPlasma;
%nodi = 200;


% l1 = 8e-2; % m (lunghezza camera 1 (quella più vicina all'upstream)
% l2 = 16e-2;
% l3 = 16e-2;
% lu = 0;
ld = 7.5e-2;
l3 = 16e-2 + ld;
lc = 1.5e-2 + l3;
l2 = 16e-2 + lc;
lb = 1.5e-2 + l2;
l1 = 8e-2 + lb;
la = 1.5e-2+l1;

%comune = 1e19;
comune0 = 1e20;
nzd = 0;
nzc = (nz2+nz3)/2;
nzb = (nz2+nz1)/2;
nza = 0;
%nz = 0;
gamma = 7;
mm = (4.002602)*1.66054*1e-27; % kg
Lpar = la; % m
%nodi = 600;
%spaceprog = 40;
spaceprog = steps;


%% costruzione funzione di radiazione
%TT_rad = logspace(-1,4,1000);
res = ne_u*tau;
%LL = F_func(TT_rad,res);

%loglog(TT_rad,LL)
linear = @(Tet,Teu) linspace(Tet,Teu,nodi);

%rad_fun = @(Tet,Teu) spline(TT_rad,LL,linear(Tet,Teu).*(linear(Tet,Teu)>=1)); % ci mette tanto
rad_fun = @(Tet,Teu) F_func(linear(Tet,Teu),res).*(linear(Tet,Teu)>=1) ; % ci mette di meno


% loglog(TT_rad,LL,'DisplayName','Dati')
% hold on
% loglog(linear(TT_rad(1),TT_rad(end)),rad_fun(TT_rad(1),TT_rad(end)),'displayName','Interpolazione')
% legend('Location','best')

%% scrittura sistema
ne_z = @(Teu,Tez) ne_u*Teu./Tez;
lunghezza = linear(0,Lpar);
nz_v = nzd*((lunghezza<ld).*(lunghezza>=0))+nz3*((lunghezza<l3).*(lunghezza>=ld))+nzc*((lunghezza<lc).*(lunghezza>=l3))+nz2*((lunghezza<l2).*(lunghezza>=lc))+nzb*((lunghezza<lb).*(lunghezza>=l2))+nz1*((lunghezza<l1).*(lunghezza>=lb))+nza*((lunghezza<la).*(lunghezza>=l1));
nzd0 = 0;
nz30 = nzPreV(3);
nz20 = nzPreV(2);
nz10 = nzPreV(1);
nzc0 = (nz20+nz30)/2;
nzb0 = (nz20+nz10)/2;
nza0 = 0;
nz_v0 = nzd0*((lunghezza<ld).*(lunghezza>=0))+nz30*((lunghezza<l3).*(lunghezza>=ld))+nzc0*((lunghezza<lc).*(lunghezza>=l3))+nz20*((lunghezza<l2).*(lunghezza>=lc))+nzb0*((lunghezza<lb).*(lunghezza>=l2))+nz10*((lunghezza<l1).*(lunghezza>=lb))+nza0*((lunghezza<la).*(lunghezza>=l1));


cz = @(Teu,Tez) nz_v./ne_z(Teu,Tez);


qqfun = @(Tet,Teu,qqt) (qqt^2+2*ke_par*ne_u^2*Teu^2*trapz(linear(Tet,Teu),rad_fun(Tet,Teu).*sqrt(linear(Tet,Teu)).*nz_v./ne_z(Teu,linear(Tet,Teu)).*(linear(Tet,Teu)>=1)))^0.5;
% qqfun = @(Tet,Teu,qqt) (qqt^2+2*ke_par*nz*ne_u^2*Teu^2*trapz(linear(Tet,Teu),rad_fun(Tet,Teu).*sqrt(linear(Tet,Teu))./ne_z(Teu,linear(Tet,Teu))))^0.5; % senza perdita sotto 1eV

qqfun_vett = @(Tet,Teu,qqt) (abs((qqt^2+2*ke_par*ne_u^2*Teu^2*trapz(linear(Tet,Teu),(tril( ones((length(linear(Tet,Teu)))) ).*linear(Tet,Teu)).^(1/2).*nz_v./ne_z(Teu,tril( ones((length(linear(Tet,Teu)))) ).*linear(Tet,Teu)).*tril( ones((length(linear(Tet,Teu)))) ).*rad_fun(Tet,Teu).*(linear(Tet,Teu)>=1),2))).^0.5)'; % restituisce un vettore di lunghezza pari a quella di linear, ogni elemento è uguale a q_par per la temperatura nel nodo i esimo
% qqfun_vett = @(Tet,Teu,qqt) (abs((qqt^2+2*ke_par*nz*ne_u^2*Teu^2*trapz(linear(Tet,Teu),(tril( ones((length(linear(Tet,Teu)))) ).*linear(Tet,Teu)).^(1/2)./ne_z(Teu,tril( ones((length(linear(Tet,Teu)))) ).*linear(Tet,Teu)).*tril( ones((length(linear(Tet,Teu)))) ).*rad_fun(Tet,Teu),2))).^0.5)'; % senza perdita sotto 1eV

% ricordati di togliere gli abs, sis principale
fun1 = @(Tet,Teu,qqt) abs((qqt^2-(qqu^2 - 2*ke_par*ne_u^2*Teu^2*trapz(linear(Tet,Teu),rad_fun(Tet,Teu).*sqrt(linear(Tet,Teu)).*nz_v./ne_z(Teu,linear(Tet,Teu)))))/qqu^2);
fun2 = @(Tet,Teu,qqt) abs((qqt - gamma*(ne_u*Teu*1.602e-19/2)*sqrt(2*Tet*1.602e-19/mm))/qqu);
fun3 = @(Tet,Teu,qqt) abs((Lpar-ke_par*trapz(linear(Tet,Teu),(linear(Tet,Teu)).^(5/2)./qqfun_vett(Tet,Teu,qqt)))/Lpar);

%altro sist boh inutile
% fun1 = @(Tet,Teu,qqt) (qqt^2-(qqu^2 - 2*ke_par*nz*ne_u^2*Teu^2*trapz(linear(max(Tet,1),Teu),rad_fun(max(Tet,1),Teu).*sqrt(linear(max(Tet,1),Teu))./ne_z(Teu,linear(max(Tet,1),Teu)))))/qqu^2;
% fun2 = @(Tet,Teu,qqt) (qqt - gamma*(ne_u*Teu*1.602e-19/2)*sqrt(2*max(Tet,1)*1.602e-19/mm))/qqu;
% fun3 = @(Tet,Teu,qqt) (Lpar-ke_par*trapz(linear(max(Tet,1),Teu),(linear(max(Tet,1),Teu)).^(5/2)./qqfun_vett(max(Tet,1),Teu,qqt)))/Lpar;

% altro altro sist (tengo q costante nella 3)
% fun1 = @(Tet,Teu,qqt) abs((qqt^2-(qqu^2 - 2*ke_par*ne_u^2*Teu^2*trapz(linear(Tet,Teu),rad_fun(Tet,Teu).*sqrt(linear(Tet,Teu)).*nz_v./ne_z(Teu,linear(Tet,Teu)))))/qqu^2);
% fun2 = @(Tet,Teu,qqt) abs((qqt - gamma*(ne_u*Teu*1.602e-19/2)*sqrt(2*Tet*(Tet>=0)*1.602e-19/mm))/qqu);
% fun3 = @(Tet,Teu,qqt) abs((Lpar-ke_par*trapz(linear(Tet,Teu),((linear(Tet,Teu).*(linear(Tet,Teu)>=0))).^(5/2))/qqfun(Tet,mean([Tet,Teu]),qqt))/Lpar);
F_sys = @(vett) [fun1(vett(1),vett(2),vett(3));fun2(vett(1),vett(2),vett(3));fun3(vett(1),vett(2),vett(3))];
%guess = [0.08,2,1e5]';
 %guess = [1.4;10;5e6]; % quello da mettere
guess = solPreV;
% guess = [9.670034050925897e-05;8.221124191752168;3.320006847362548e+04];
%% fsolve

 %SOLUZIONE = (fsolve(F_sys,guess))
 %SOLUZIONE = lsqnonlin(F_sys,guess)

%  progr = linspace(comune0,max(nz_v),10);
prog = zeros(spaceprog,nodi);
% for ii = 1:length(nz_v)
%     prog(:,ii) = logspace(log10(nz_v0(ii)),log10(nz_v(ii)),spaceprog);
% end
guessM = zeros(3,spaceprog);
if guess(1) >=1

    for jj = 1:spaceprog
        
        %nz_v = prog(jj,:);
    nz_v = nzd*((lunghezza<ld).*(lunghezza>=0))+nz3*((lunghezza<l3).*(lunghezza>=ld))+nzc*((lunghezza<lc).*(lunghezza>=l3))+nz2*((lunghezza<l2).*(lunghezza>=lc))+nzb*((lunghezza<lb).*(lunghezza>=l2))+nz1*((lunghezza<l1).*(lunghezza>=lb))+nza*((lunghezza<la).*(lunghezza>=l1));
    cz = @(Teu,Tez) nz_v./ne_z(Teu,Tez);
    
    
    qqfun = @(Tet,Teu,qqt) (qqt^2+2*ke_par*ne_u^2*Teu^2*trapz(linear(Tet,Teu),rad_fun(Tet,Teu).*sqrt(linear(Tet,Teu)).*nz_v./ne_z(Teu,linear(Tet,Teu)).*(linear(Tet,Teu)>=1)))^0.5;
    % qqfun = @(Tet,Teu,qqt) (qqt^2+2*ke_par*nz*ne_u^2*Teu^2*trapz(linear(Tet,Teu),rad_fun(Tet,Teu).*sqrt(linear(Tet,Teu))./ne_z(Teu,linear(Tet,Teu))))^0.5; % senza perdita sotto 1eV
    
    qqfun_vett = @(Tet,Teu,qqt) (abs((qqt^2+2*ke_par*ne_u^2*Teu^2*trapz(linear(Tet,Teu),(tril( ones((length(linear(Tet,Teu)))) ).*linear(Tet,Teu)).^(1/2).*nz_v./ne_z(Teu,tril( ones((length(linear(Tet,Teu)))) ).*linear(Tet,Teu)).*tril( ones((length(linear(Tet,Teu)))) ).*rad_fun(Tet,Teu).*(linear(Tet,Teu)>=1),2))).^0.5)'; % restituisce un vettore di lunghezza pari a quella di linear, ogni elemento è uguale a q_par per la temperatura nel nodo i esimo
    % qqfun_vett = @(Tet,Teu,qqt) (abs((qqt^2+2*ke_par*nz*ne_u^2*Teu^2*trapz(linear(Tet,Teu),(tril( ones((length(linear(Tet,Teu)))) ).*linear(Tet,Teu)).^(1/2)./ne_z(Teu,tril( ones((length(linear(Tet,Teu)))) ).*linear(Tet,Teu)).*tril( ones((length(linear(Tet,Teu)))) ).*rad_fun(Tet,Teu),2))).^0.5)'; % senza perdita sotto 1eV
    
        fun1 = @(Tet,Teu,qqt) abs((qqt^2-(qqu^2 - 2*ke_par*ne_u^2*Teu^2*trapz(linear(Tet,Teu),rad_fun(Tet,Teu).*sqrt(linear(Tet,Teu)).*nz_v./ne_z(Teu,linear(Tet,Teu)))))/qqu^2);
    fun2 = @(Tet,Teu,qqt) abs((qqt - gamma*(ne_u*Teu*1.602e-19/2)*sqrt(2*Tet*1.602e-19/mm))/qqu);
    fun3 = @(Tet,Teu,qqt) abs((Lpar-ke_par*trapz(linear(Tet,Teu),(linear(Tet,Teu)).^(5/2)./qqfun_vett(Tet,Teu,qqt)))/Lpar);
    F_sys = @(vett) [fun1(vett(1),vett(2),vett(3));fun2(vett(1),vett(2),vett(3));fun3(vett(1),vett(2),vett(3))];
    
    
    
        %% newton
        jfun2a=@(x) numerical_jacobian(F_sys,x,1e-5);
        if jj == length(1:spaceprog)
            toll = 1e-2;
            max_iter = 500;
            
        else
            toll=1e-1;
            max_iter = 50;
            
        end
        %max_iter=500;
        nz_v;
        if flagShow == 1
      jj
        end
    %   guess(1)
    %   guess(2)
    %   guess(3)
        
            [PROVA,err,residual,niter]=myNewton_Jac(F_sys,jfun2a,guess,toll,max_iter);
            flag_det = 0;
        
            
    
              
    
        
    %   if jj > 2
    %       guess = abs(guessM(:,jj-2)+(guessM(:,jj-1)-guessM(:,jj-2))/(max(prog(jj-1,:))-max(prog(jj-2,:)))*(max(prog(jj,:)-max(prog(jj-2,:)))))
    %   else
      guess = PROVA; 
    %   end
      guessM(:,jj) = guess;
      
    end

else

% % sistema detachment
%     qqd = guess(3);
%     flag_det = 1;
%     nz_aux = @(lunghezza) nzd*((lunghezza<ld).*(lunghezza>=0))+nz3*((lunghezza<l3).*(lunghezza>=ld))+nzc*((lunghezza<lc).*(lunghezza>=l3))+nz2*((lunghezza<l2).*(lunghezza>=lc))+nzb*((lunghezza<lb).*(lunghezza>=l2))+nz1*((lunghezza<l1).*(lunghezza>=lb))+nza*((lunghezza<la).*(lunghezza>=l1));
%     nz_z = @(zz) nz_aux(zz);
%     nz_Te = @(Tez,Teu,Ld) nz_z( (Tez-Teu)*Ld/(1-Teu) );
%     
%     qqfund = @(Teu,Tez,Ld) (qqd^2+2*ke_par*ne_u*Teu*trapz(linear(1,Tez),rad_fun(1,Tez).*nz_Te(linear(1,Tez),Teu,Ld).*(linear(1,Tez)).^(3/2).*(linear(1,Tez)>=1)))^0.5;
%     
%     qqfun_vettd = @(Teu,Tez,Ld) (abs((qqd^2+2*ke_par*ne_u^2*Teu^2*trapz(linear(1,Tez),(tril( ones((length(linear(1,Tez)))) ).*linear(1,Tez)).^(1/2).*nz_Te(linear(1,Tez),Teu,Ld)./ne_z(Teu,tril( ones((length(linear(1,Tez)))) ).*linear(1,Tez)).*tril( ones((length(linear(1,Tez)))) ).*rad_fun(1,Tez).*(linear(1,Tez)>=1),2))).^0.5)'; % restituisce un vettore di lunghezza pari a quella di linear, ogni elemento è uguale a q_par per la temperatura nel nodo i esimo
%     
%     fun1d = @(Teu,Ld) abs(qqu^2 - qqfund(Teu,Teu,Ld)^2)/qqu^2;
%     fun2d = @(Teu,Ld) abs((Ld-ke_par*trapz(linear(1,Teu),(linear(1,Teu)).^(5/2)./qqfun_vettd(Teu,Teu,Ld))))/Lpar;
%     
%     
%     F_sysd = @(vett) [fun1d(vett(1),vett(2));fun2d(vett(1),vett(2))];
% 
%     SOL = fsolve(F_sysd,guess);  % devi sistemare il guess, devi mettere [guess(2);Lpar] però poi deve cambiare ogni volta perchè devi cambiare proprio il solPreV di wall
%     PROVA = [1;SOL(1);SOL(2)];
% %% distaccamento
% fun4 = @(Teu) qqu^2-2*ke_par*ne_u*Teu*trapz(linear(0,Teu),rad_fun(0,Teu).*(linear(0,Teu)).^(5/2).*nz_v);
% 
% Teu_dist = fzero(fun4,9);
% 
% intgrRad = qqfun_vett(0,Teu_dist,0.1);
% 
% %Ldist = ke_par*trapz(linspace(1,Teu_dist,length(intgrRad(find(intgrRad):end))),(linspace(1,Teu_dist,length(intgrRad(find(intgrRad):end)))).^(5/2)./(intgrRad(find(intgrRad):end))) 
% Ldist = ke_par*trapz(linspace(0,Teu_dist,length(intgrRad(find(intgrRad):end))),(linspace(1,Teu_dist,length(intgrRad(find(intgrRad):end)))).^(5/2)./(intgrRad(find(intgrRad):end))); % senza perdita sotto 1eV

end
%% calcolo potenze
PirrTot = Pot*(qqu-PROVA(end))/qqu;

Tet = PROVA(1);
Teu = PROVA(2);
qqt = PROVA(3);



T_vett = linear(Tet,Teu);
lungh_d_v = lunghezza<=ld;
ii_d = sum(lungh_d_v);
% ii_3 = sum((lunghezza<=l3).*(lunghezza>ld));
% ii_2 = sum((lunghezza<=l2).*(lunghezza>l3));
% ii_1 = sum((lunghezza<=l1).*(lunghezza>l2));
ii_3 = sum(lunghezza<=l3);
ii_2 = sum(lunghezza<=l2);
ii_1 = sum(lunghezza<=l1);

% qq_d = qqfun(T_vett(1),T_vett(ii_d),qqt);
% qq_3 = qqfun(T_vett(1),T_vett(ii_3),qqt);
% qq_2 = qqfun(T_vett(1),T_vett(ii_2),qqt);
% qq_1 = qqfun(T_vett(1),T_vett(ii_1),qqt);

qqv = qqfun_vett(Tet,Teu,qqt);
qq_d = qqv(ii_d);
qq_3 = qqv(ii_3);
qq_2 = qqv(ii_2);
qq_1 = qqv(ii_1);

Pirr1 = Pot*(qq_1-qq_2)/qq_1;
Pirr2 = (Pot-Pirr1)*(qq_2-qq_3)/qq_2;
Pirr3 = (Pot-Pirr1-Pirr2)*(qq_3-qq_d)/qq_3;

potenze_vettore = [Pirr1,Pirr2,Pirr3]';
soluzione = PROVA;
flagMaxIter = 0;
if niter == max_iter
    flagMaxIter = 1;
end
AA = [potenze_vettore,soluzione,[flagMaxIter;niter;max_iter] ];


%plot(linear(Tet,Teu),qqfun_vett(Tet,Teu,qqt))






return
% %% punto fisso
% 
% iter = 1;
% toll = 1e-4;
% err = toll+1;
% iterMax = 200;
% solOld = guess;
% solNew = solOld;
% while iter < iterMax && err > toll 
% %     solNew(2) = (qqu^2-solOld(3)^2)/(2*ke_par*nz*ne_u*trapz(linear(solOld(1),solOld(2)),rad_fun(solOld(1),solOld(2)).*(linear(solOld(1),solOld(2))).^(3/2) ));    
% %     solNew(3) = gamma*ne_u*(solOld(2)*1.602e-19)/2*sqrt(2*(solOld(1)*1.602e-19)/mm);
% %     fun3_aux = @(Tet) -Lpar+ke_par*trapz(linear(abs(Tet),solOld(2)),(linear(abs(Tet),solOld(2))>=1).*(linear(abs(Tet),solOld(2))).^(3/2)./qqfun_vett(abs(Tet),solOld(2),solOld(3)));
% %     solNew(1) = fzero(fun3_aux,solOld(1));
% 
%     %provo altro sis
%     vv = linear(solOld(1),solOld(2));
%     ss = linear(0,Lpar);
%     ww = find(vv>=1);
%     
%     %cose strane
%     %fun4_aux = @(Teu) -Teu + vv(end-1)+2*(qqfun(solOld(1),solOld(2),solOld(3))/ke_par*(ss(end)-ss(end-1))/(vv(end)-vv(end-1)))^(2/5);
%     %solNew(2) = vv(end-1)+2*(qqfun(solOld(1),solOld(2),solOld(3))/ke_par*(ss(end)-ss(end-1))/(vv(end)-vv(end-1)))^(2/5);
% 
%     %solNew(2) = (qqfun(solOld(1),solOld(2),solOld(3))/ke_par*(ss(end)-ss(end-1))/(vv(end)-vv(end-1)))^(2/5); % metti questa
%     solNew(2) = (qqfun(solOld(1),solOld(2),solOld(3))/ke_par*(ss(end)-ss(end-1))/(vv(end)-vv(end-1)))^(2/5);
%     solNew(1) = (solOld(3)/(gamma*ne_u*(solOld(2)*1.602e-19)/2))^2*mm/2/1.602e-19; %inizio
%     if 2*ke_par*nz*ne_u*solOld(2)*trapz(linear(solOld(1),solOld(2)),rad_fun(solOld(1),solOld(2)).*(linear(solOld(1),solOld(2)).*(linear(solOld(1),solOld(2))>=1)).^(3/2) ) >qqu^2
%         solNew(3) = 0;
%     else
%         solNew(3) = (qqu^2-2*ke_par*nz*ne_u*solOld(2)*trapz(linear(solOld(1),solOld(2)),rad_fun(solOld(1),solOld(2)).*(linear(solOld(1),solOld(2)).*(linear(solOld(1),solOld(2))>=1)).^(3/2) ))^0.5; %nell'integr
%     end
%     err = norm(solNew-solOld)/norm(solNew);
%     iter = iter+1;
%     solOld = solNew;
% end
% 
% Lreal = Lpar-ss(ww(1));
% fun3Real = @(Tet,Teu,qqt) (Lreal-ke_par*trapz(linear(Tet,Teu),(linear(Tet,Teu)).^(5/2)./qqfun_vett(Tet,Teu,qqt)))/Lreal;






