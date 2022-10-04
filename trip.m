%function [Fobj] = trip(va) % comment for direct use

% Validations
%               c_trip=0    V=32.4km/h d=1000m p=0% (m=78kg, Cr=0.0059, SCx=0.43*0.8) => P=194W E=6Wh Paero=79% OK
%               c_trip=1    V=25.0km/h d=1000m p=2% (m=78kg, Cr=0.0059, SCx=0.43*0.8) => P=211W E=8Wh Paero=34% OK
%               c_trip=2    Accélération Vmax=52km/h in 200m p=0% (Cr=0.0059, SCx=0.43*0.8) =>                  NO
                            % mv=12 => P=917W E=4Wh Paero=34% Pi=60%
                            % mv= 8 => P=888W E=4Wh Paero=36% Pi=59%
%               c_trip=6    Ascension du Ventoux (78,0.005,0.43*0.80,0,0) P=Ph+Pe=200+0 xzP=1 t=1h55mins V=11.8km/h E=380Wh
%               c_trip=9    500km 50km/h 10h VCCAES v0 plat à V=Cte / t1.m E = 5000Wh or 500W x 10h
 
% Compute VCCAES record 500km in a day with Solar Panels
% 500km D+1500m 10h 50km/h  with VCCAES     v0 Eb0=1400Wh PV=1.6m2 (m,Cr,SCx)=(110,0.0080,0.26) le 20 Juin d=400km Eb(end)= -418Wh
%                           with Solarboost v0 Eb0=3000Wh PV=2.0m2 (m,Cr,SCx)=(110,0.0080,0.32) le 20 Juin d=500km Eb(end)=  597Wh

clear all
addpath(genpath(cd));

rho      = 1.2;
g        = 9.8;
eta_m    = 0.90;    % engine  efficiency
eta_mppt = 0.97;    % mppt    efficiency
eta_d    = 0.95;    % driving efficiency
tp       = 10;      % Pilote behavior / acceleration (start of each segment)

% Env. data (Sun & Wind conditions)
lat     = 43.57; long = -1.47;                          % Location : ISAE-SUPAERO
Jour    = 20; Mois = 6; Annee=2023; t_start = 7;        % Date : Départ à 7h UTC le 20 Juin 2023 (voir coeff)
v_wind  = 0/3.6;                                        % wind speed in the opposite direction

% Trip choice
%c_trip  = 6; mode=2;xzP=0;zg=1;    % Trip 5,6,7        = (x,z,v)       Ascension Ventoux à V fixée avec moteur Pm(t).
%c_trip  = 6; mode=2;xzP=1;zg=1;    % Trip 5,6,7        = (x,z,ph,pe)   Ascension Ventoux à P=Ph+Pe=Cte.
%c_trip  = 8; mode=2;xzP=1;zg=1;    % Trip 8            = gpx file
 c_trip  = 9; mode=2;xzP=0;zg=1;    % Trip 9,10,11, ... = (x,z,v)       VCCAES    500km

%% Choose Bike + Assistance mode for VAE
%       mv=  14; Cr=0.0070; mc  = 70; SCx = 0.56*0.80;          S_PV=0.0; Eb0=   220; Pht=200; Pet=0.66*Pht;    % VAE D31
%       mv=  80; Cr=0.0090; mc  = 70; SCx = 0.85*1.00*1.0*0.6;  S_PV=2.0; Eb0=  3000; Pht=200; % VCCAES SolarBoost v00
        mv=  40; Cr=0.0080; mc  = 70; SCx = 0.80*1.00*1.0*0.4;  S_PV=2.0; Eb0=  3000; Pht=200; % VCCAES SolarBoost v0
%       mv=  30; Cr=0.0060; mc  = 70; SCx = 0.80*0.80*1.0*0.3;  S_PV=2.0; Eb0=  3000; Pht=100; % VCCAES SolarBoost v1
Eb(1)   = Eb0;
m       = mv+mc;

%% Data
switch c_trip
    case 0 % Validation V=32.4km/h d=1000m p=0% (m=78kg, Cr=0.0059, SCx=0.43*0.8) => P=194W E=6Wh Paero=79%
        t=[0 110];z=[0 0];v=[9 9];
    case 1 % Validation V=25.0km/h d=1000m p=2% (m=78kg, Cr=0.0059, SCx=0.43*0.8) => P=211W E=8Wh Paero=34%
        t=[0 143]; z=[0 20]; v=[7 7];
    case 2 % Validation accélération Vmax=52km/h in 200m p=0% (Cr=0.0059, SCx=0.43*0.8) => 
        % mv=12 => P=917W E=4Wh Paero=34% Pi=60%
        % mv= 8 => P=888W E=4Wh Paero=36% Pi=59%
        t=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];z=zeros(1,length(t));v=[0 3 5 6.7 7.9 8.8 9.5 10.2 10.9 11.5 12.1 12.5 12.9 13.25 13.5 13.85 14.15 14.4];
    case 3 % Démarrage feu rouge V=0 à 28km/h en 10s soit a=0.78m/s2 p=0% => P=379W Paero=13% Ei=80%
        t=[0 2 4 6 8 10];z=[0 0 0 0 0 0];v=[0 3.0 4.7 6.0 7.0 7.8];
    case 4 % Démarrage feu rouge V=0 à 36km/h en 10s soit a=1.00m/s2 p=0% => P=450W Paero=17% Ei=78%
        t=[0 2 4 6 8 10];z=[0 0 0 0 0 0];v=[0 3.8 6.0 7.7 9.0 10];
    case 5 % Montée 500m à 2% et descente 500m à -2% => P=190W Paero=83%
        x=[0 500 1000];z=[0 10 0];v=[10 6.9 10.5];
    case 6 % Ascension Mont Ventoux
        cv=1;
        x=[  0 400 800 1800 2800 3800 4800 5800 6800 7800 8800 9800 10800 11800 12800 13800 14800 15800 16800 17800 18800 19800 21800 22520];   % [m]   GPS
        z=[326 320 311  335  371  415  473  529  581  690  785  890   982  1068  1154  1236  1315  1396  1461  1515  1591  1662  1837  1893];   % [m]   altitude GPS
        v=[  5   8  11   10    9    8  7.5    7  6.5  5.0  5.0  5.0     5     5     5     5     4     4     4     4     4     4     4     5]*cv;% [m/s] vitesse  GPS
    case 8 % [x z v] = f(gpx file) created by pGPX.m
         load('parcours.mat')           % load [xzv]
        v=v/3.6; % from km/h to m/s
    case 9 % 500km D+0m 10h - Validation / t1.m OK
        vm = 13.9;
        x=[  0 100 5000  20000 30000 40000 50000 60000 70000 80000 120000 150000 200000 250000 300000 350000 400000 450000 480000 490000  500000];   % [m]   GPS
        z=[  0   0    0      0     0     0     0     0     0     0      0      0      0      0      0      0      0      0      0      0       0];   % [m]   altitude GPS
        v=[  0   7   vm     vm    vm    vm    vm    vm    vm    vm     vm     vm     vm     vm     vm     vm     vm     vm     vm     vm      vm];   % [m/s] vitesse  GPS
end

%% From [xzv] to [xzv] with smoothing and few points
%v = smoothdata(v,'gaussian',20); % to discard creazy points of GPS

%% [x,z,v] = f([x,z], Pht+Pet)  (t.m)
if xzP==1 % v=f(P)
    fprintf(' from [xz] or [xzv] to [xzP]'); 
    for i = 2 : length(z) % loop on segments
    %     Ph_target = 200;Pm_target = 0.8*Pmax_engine;coeff = 0.8;Ps = 1000*0.22*S_PV.*coeff; % Flux solaire récupéré [W]
        v0        = 7;
        P_target  = Pht+Pet;
        a         = 0.01; % arbitrary acceleration
        p         = (z(i)-z(i-1))/(x(i)-x(i-1)); % slope
        f = @(x,rho,SCx,Cr,p,m,g,a,eta_d,P_target)(0.5*rho*x^3*SCx+(Cr+p)*m*g*x+m*a*x)/eta_d-P_target; % v = f(P_target, ...)
        options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter');
        fun = @(x)f(x,rho,SCx,Cr,p,m,g,a,eta_d,P_target);
        [v(i),fval,exitflag,output] = fsolve     (fun,v0,options)   % x / f(x)=0
    end
end

%% Init
x (1) = 0;
t (1) = 0;dt(1)=0;
a (1) = 0;
vn(1) = v(1);
dz    = 0;

%% Eq of mvt a(t)=f(v,t)
switch mode
    case 1 % t,z,v given : compute (x,a,p)=f(t,z,v)
        for i = 2 : length(z)                          % for each segment
            dt(i)    = t(i)-t(i-1);                    % ?
            x (i)    = x (i-1) + v(i)*dt(i);
            a (i)    = (v(i)-v(i-1)) / dt(i);          % v = ve = v - va, P(v)
            vn(i)    = vn(i-1) + a(i)*dt(i);           % for validation vn=v
            p (i)    = (z(i)-z(i-1)) / (x(i)-x(i-1));  % slope in %
            dz       = dz + max(z(i)-z(i-1),0);        % compute D+
            dts      = min(tp,dt(i)/10);               % time for acceleration
            as(i)    = (vn(i)-vn(i-1))/dts;            % initial acc (pilote behavior)
            coeff(i) = 0;                              % no sun, no PV
        end
        fprintf(' mean dt = %3.1f sec \n',mean(dt));
    case 2 % x,z,v given : compute (t,a,p) = f(x,z,v) - ebikemotion D31
        fprintf(' mode 2 : compute (t,a,p) = f(x,z,v) - ebikemotion, Komoot \n');
        for i = 2 : length(z)                       % for each segment
            t (i) = t (i-1) + (x(i)-x(i-1))/v(i);   % t=f(x,v) if a=cte
            dt(i) = t(i)-t(i-1);                    % time
            a (i) = (v(i)-v(i-1)) / dt(i);          % a=f(v,t) mean acceleration
            vn(i) = vn(i-1) + a(i)*dt(i);           % for validation vn=v
            p (i) = (z(i)-z(i-1)) / (x(i)-x(i-1));  % slope in %
            dz    = dz + max(z(i)-z(i-1),0);        % compute D+            
            dts   = min(tp,dt(i)/2);                % acceleration on half of the segment
            as(i) = (vn(i)-vn(i-1))/dts;            % initial acc (pilote behavior)
            % compute solar coeff & power at location (lat, long) at timet / PV 0°.
            PV_angle(i)  = 0; % Angle between PV and Ground
            mDateVec  = datenum([Annee,Mois,Jour,t_start+t(i)/3600,0,0]); % Heure UTC
            [wAz,wEl] = SolarAzEl(mDateVec,zeros(size(mDateVec,1),1)+lat,zeros(size(mDateVec,1),1)+long,zeros(size(mDateVec,1),1));
            if wEl < 0
                ang = 0;    % sun below horizon
            else
                ang = wEl+PV_angle(i);
            end
            coeff(i)  = cos((90-ang)*pi/180);coeff=max(coeff,0); % cos(theta_loc) : 
        end
end        

%% Forces
Fa = 0.5*1.225.*v.^2*SCx;   % Aero     force
Fr = Cr*m*9.8.*cos(p);      % Rolling  force
Fg =    m*9.8.*sin(p);      % Gravity  force
Fi =    m.*as;              % Inertial force (start/pilote behavior)
%Fi =   m.*a;               % Inertial force (mean)
F  = Fa + Fr + Fg + Fi;     % Total force

Fa(1)=0;Fr(1)=0;Fg(1)=0;Fi(1)=0; % Correct initial point

%% Power consumption
Pa = Fa.*(v+v_wind);        % Aero power
Pr = Fr.*(v+v_wind);        % Rolling power
Pg = Fg.*(v+v_wind);        % Gravity power
Pi = Fi.*(v+v_wind);        % Inertial power

%% Power received
Ps = 1000*0.22*S_PV.*coeff; % Flux solaire récupéré [W] = 1000W/m2 * rendement PV * surface * coeff(place,date,hour)

switch xzP
    case 0 % Ph = Pht
        Ph = ones(1,length(t))*Pht; % Human power Ph=Pht only true if v not given
        %Ph = F.*(v+v_wind); % Ph+Pe = F.V
    case 1 % Ph = Pht
        Ph = ones(1,length(t))*Pht; % Human power Ph=Pht only true if v not given
end

Pc = Pa + Pr + Pg + Pi;     % Power needed
Pp = Ph + Ps;               % Power received (human + solar)
Pm = max(0,Pc - Pp);        % Power consumption electric = Motor power
%Pm = max(0,Pc - Ph);       % Power consumption electric = Motor power

Pa(1)=0;
Pr(1)=0;
Pg(1)=0;
Pi(1)=0;
Ph(1)=0;
Pm(1)=0;
Pc(1)=0;
Ps(1)=0; 

PC(:,1)=Pa(:);PC(:,2)=Pr(:);PC(:,3)=Pg(:);PC(:,4)=Pi(:);PP(:,1)=Ph(:);PP(:,2)=Pp(:);  % Bar plot power production

%% Cummulated Energy [Wh]
E=0;Ea=0;Er=0;Eg=0;Ei=0;Eh=0;Es=0;Em=0;El=0;
for i = 2 : length(t)
    %dt    = t(i)-t(i-1);                        % segment dt
    Ec(i) = (Pc(i)/eta_d)*dt(i)/3600;            % Full     energy consumption = (Ph+Pe) *dt
    Ea    = Ea + Pa(i)*dt(i)/3600;               % Aero     energy consumption
    Er    = Er + Pr(i)*dt(i)/3600;               % Rolling  energy consumption
    Eg    = Eg + max(Pg(i),0)*dt(i)/3600;        % Gravity  energy consumption (no gravity energy recup)
    Ei    = Ei + max(Pi(i),0)*dt(i)/3600;        % Inertial energy consumption
    Eh    = Eh + Ph(i)*dt(i)/3600;               % Human    energy received
    Es    = Es + Ps(i)*dt(i)/3600;               % Solar    energy received
    Em    = Em + Pm(i)*dt(i)/3600;               % Electric energy consumption
    Ef    = (-Pm(i)/(eta_m*eta_mppt))*dt(i)/3600;       % Energy flux out of battery
%   Ef    = (-Pm(i)/eta_m + Ps(i)*eta_mppt)*dt(i)/3600; % Energy flux out of battery
    Eb(i) = min( Eb0 , Eb(i-1)+Ef );            % Battery energy stock
    Elm   = (Pm(i)/eta_m - Pm(i) )*dt(i)/3600;  % Energy losses / motor efficiency
    Eld   = (Pc(i)/eta_d - Pc(i) )*dt(i)/3600;  % Energy losses / mechanical efficiency
    Els   = (Ps(i)-eta_mppt*Ps(i))*dt(i)/3600;  % Energy losses / mppt efficiency
    El    = El + ( Elm + Els + Eld );           % Cumulative Energy losses
    E     = E  + Ec(i);                         % Total energy consumption
end

%% Plot
if zg==1
    
    figure(1);
    nl = 6;
    subplot(nl,1,1);plot(x/1000,z,'o-')    ;ylabel('z[m]');axis([0 x(end)/1000 min(z)-1 max(z)+1]);
    subplot(nl,1,2);plot(x/1000,p*100,'o-');ylabel('slope[%]');axis([0 x(end)/1000 min(p*100)-1 max(p*100)+1]);
    subplot(nl,1,3);plot(x/1000,v*3.6,'o-');ylabel('Speed [km/h]');axis([0 x(end)/1000 min(v*3.6) max(v*3.6)]);
    %subplot(nl,1,3);plot(x/1000,Pm,'o-');ylabel('Engine power [W]');
    subplot(nl,1,4);plot(x/1000,a,'o-'    );ylabel('a[m/s2]');axis([0 x(end)/1000 min(a) max(a)]);
    %subplot(nl,1,4);plot(t,coeff,'o-'    );xlabel('t[s]');ylabel('Solar coeff');
    subplot(nl,1,5);plot(x/1000,Pc,'o-',x/1000,Pp,'o-',x/1000,Ps,'o-');ylabel('Power [W]');legend('Pc','Pp=Ph+Ps','Ps');
    subplot(nl,1,6);plot(x/1000,Eb,'o-');xlabel('x[km]');ylabel('Eb[Wh]'); axis([0 x(end)/1000 0 Eb(1)+0.01]);
    %hold on;

    figure(2);
    nl = 10;
    subplot(nl,1,1);bar(z)    ;ylabel('z [m]');axis([-inf inf 0 inf]);
    %text(0,1500,num2str(round(p*100,1)));
    subplot(nl,1,2);bar(p*100) ;ylabel('slope[%]');
    subplot(nl,1,3);bar(v*3.6) ;ylabel('V [km/h]');axis([-inf inf 0 inf]);
    subplot(nl,1,4);bar(dt/60) ;ylabel('t [min]');axis([-inf inf 0 inf]);
    %subplot(nl,1,4);bar(a)     ;ylabel('Acc.  [m/s2]');
    subplot(nl,1,5);bar(Ph)    ;ylabel('Ph [W]');axis([-inf inf 0 inf]);
    subplot(nl,1,6);bar(Ps)    ;ylabel('Ps [W]');axis([-inf inf 0 inf]);
    subplot(nl,1,7);bar(Pm)    ;ylabel('Pm [W]');axis([-inf inf 0 inf]);
    subplot(nl,1,8);bar(PC,'stacked');ylabel('Psink   [W]');legend('Pa','Pr','Pg','Pi');grid on;axis([-inf inf 0 inf]);
    %subplot(nl,1,4);bar(PP,'stacked');ylabel('Psource [W]');legend('Ph','Pm');grid on;axis([-inf inf 0 inf]);
    %subplot(nl,1,6);bar(Pp-Pc,'stacked');ylabel('Pflux[W]');legend('Power flux');grid on;axis([-inf inf -inf inf]);
    subplot(nl,1,9);bar(Ec,'stacked');xlabel('segments');ylabel('E[Wh]');axis([-inf inf 0 inf]);
    %subplot(nl,1,6);bar(Pc.*dt/3600,'stacked');xlabel('segments');ylabel('Energy[Wh]');
    subplot(nl,1,10);bar(Eb);ylabel('Eb[Wh]');axis([-inf inf 0 inf]);

    % figure(3);
    % nl = 6;
    % subplot(nl,1,1);bar(p*100);xlabel('segments');ylabel('Slope [%]');
    % subplot(nl,1,2);bar(dt)   ;xlabel('segments');ylabel('Time [s]');
    % subplot(nl,1,3);bar(v*3.6);xlabel('segments');ylabel('Speed [km/h]');
    % subplot(nl,1,4);bar(a)    ;xlabel('segments');ylabel('Acc.  [m/s2]');
    % subplot(nl,1,5);b=bar(P)  ;xlabel('segments');ylabel('Power [W]');
    % subplot(nl,1,6);bar(PP)   ;xlabel('segments');ylabel('Power [W]');legend('Pa','Pr','Pg','Pi');
end

%% Range estimate
ir=0;
for i = 1 : length(Eb)-1
    if Eb(i)>0 && Eb(i+1)<0
        ir = i;
    end
end
if ir ==0
    ir=length(Eb);
end
Range = x(ir);
% Duration estimate
if t(end)/60 > 60
    Hours   = floor(t(end)/3600);
    Minutes = t(end)/60-Hours*60;
else
    Hours   = 0;
    Minutes = t(end)/60;
end

%% Print
d = x(end)/1000;
fprintf(' Lat %4.1fN Long %4.1fE  Date %i/%i/%i  Start at %ih \n',lat,long,Jour,Mois,Annee,t_start);
fprintf(' Bike (m,Cr,SCd,Eb,S_pv)=(%3.0f,%5.4f,%4.3f,%4.0f,%4.1f) \n',mv+mc,Cr,SCx,Eb0,S_PV);
fprintf(' Trip %5.0fkm  D+%5.0fm %3.0fh%2.0fmins %4.1fkm/h  Eh=%4.0fWh \n',x(end)/1000,dz,Hours,Minutes,x(end)/t(end)*3.6,Eh);
fprintf('-------------------------------------------------------------\n');
fprintf('            min    moy    max \n');
fprintf(' Speed   %6.0f %6.0f %6.0f km/h  \n',min (v(2:end)*3.6),mean(v(2:end)*3.6),max (v(2:end)*3.6));
fprintf(' Pc      %6.0f %6.0f %6.0f W     \n',min (Pc(2:end)),mean(Pc(2:end)),max (Pc(2:end))); % to generalize...
fprintf(' Ph      %6.0f %6.0f %6.0f W     \n',min (Ph(2:end)),mean(Ph(2:end)),max (Ph(2:end))); % to generalize...
fprintf(' Ps      %6.0f %6.0f %6.0f W     \n',min (Ps(2:end)),mean(Ps(2:end)),max (Ps(2:end))); % to generalize...
fprintf(' Pm      %6.0f %6.0f %6.0f W     \n',min (Pm(2:end)),mean(Pm(2:end)),max (Pm(2:end))); % to generalize...
fprintf('           %3.0f%%   %3.0f%%   %3.0f%% Aerodynamic power  \n',min(Pa(2:end))/min(Pc(2:end))*100,mean(Pa(2:end))/mean(Pc(2:end))*100,max(Pa(2:end))/max(Pc(2:end))*100);
fprintf('-------------------------------------------------------------\n');
fprintf(' Energy E/d = Eh+Es+Eb-El = %4.0f = %4.0f + %4.0f + %4.0f - %4.0f Wh/km \n',E/d,Eh/d,Es/d,(Eb0-Eb(end))/d,El/d);
fprintf(' Energy received    = Eh+Es+Eb-El = %4.0f = %4.0f + %4.0f + %4.0f - %4.0f Wh \n'   ,E  ,Eh  ,Es  , Eb0-Eb(end)   ,El);
fprintf(' Energy consumption = Ea+Er+Eg+Ei = %4.0f = %4.0f + %4.0f + %4.0f + %4.0f Wh \n'   ,Ea+Er+Eg+Ei,Ea,Er,Eg,Ei);
fprintf(' Energy received    = Eh+Es+Eb-El = %3.0f%% = %3.0f%% +  %2.0f%% +  %2.0f%% -  %2.0f%% \n',100,100*Eh/E,100*Es/E,100*(Eb0-Eb(end))/E,100*El/E);
fprintf(' Energy consumption = Ea+Er+Eg+Ei = %3.0f%% = %3.0f%% +  %2.0f%% +  %2.0f%% +  %2.0f%% \n',100,100*Ea/E,100*Er/E,100*Eg/E,100*Ei/E);
fprintf(' Battery energy at trip end =%4.0fWh, %2.0f%% of its capacity \n',Eb(end),100*Eb(end)/Eb0);
%fprintf(' Range %4.0fkm %2.0f%% E_Flux/Stock Battery Eb(end)=%4.0fWh or %2.0f%% \n',Range/1000,100*Es/Em,Eb(end),100*Eb(end)/Eb0);

save trip.mat

