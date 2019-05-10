% Thomas R. Shannon, Fei Wang, Jose Puglisi, Christopher Weber, Donald M. Bers
% "A Mathematical Treatment of Integrated Ca Dynamics Within the Ventricular Myocyte",
% Biophys J. 2004 Nov;87(5):3351-71.

% Note:
% correction in Ito as described in Erratum in Biophys J. 2012 Apr 18;102(8):1996-2001
% correction of J_CaB_cytosol (4 Aug 2015, SM)

function output = rabbit()

clear all
close all
clc

%% Initial conditions
%load yfinal_rabbit_05      % 0.5 Hz
load yfinal_rabbit_1        % 1 Hz
%load yfinal_rabbit_2       % 2 Hz

y0 = yfinal;

%% Single Run Simulation
tspan = [0; 3e3]; % duration (ms)
HR = 1; % stimulation frequency (Hz)
options = odeset('RelTol',1e-5,'MaxStep',1,'Stats','on'); 
p = [HR];  % Parameter array for passing nondefault conditions
[t,y] = ode15s(@f,tspan,y0,options,p);

%% Final conditions
yfinal = y(end,:);
output = yfinal;

%save yfinal_rabbit_05 yfinal
%save yfinal_rabbit_1 yfinal
%save yfinal_rabbit_2 yfinal

%% Plot results
figure(1);
subplot(3,1,1); plot(t,y(:,39)); ylabel('Vm (mV)');
subplot(3,1,2); plot(t,y(:,38)); ylabel('[Ca]i (mM)');
subplot(3,1,3); plot(t,y(:,34)); ylabel('[Na]i (mM)');
xlabel('Time (ms)')

currents = calcCurrents(t,y,p);

figure(2);
subplot(3,1,1); plot(t,y(:,39)); ylabel('Vm (mV)');
subplot(3,1,2); plot(t,currents(:,1)); ylabel('ICa tot (A/F)')
subplot(3,1,3); plot(t,currents(:,2)); ylabel('INa tot (A/F)')
xlabel('Time (ms)')

function output = f(t,y,p,runType)

%% State variables
% y(1) INa - gate m
% y(2) INa - gate h
% y(3) INa - gate j 
% y(4) ICa - gate d
% y(5) ICa - gate f
% y(6) ICa - fCa_junc
% y(7) ICa - fCa_sl
% y(8) Itos - gate x
% y(9) Itos - gate y
% y(10) Itof - gate x
% y(11) Itof - gate y
% y(12) IKr - gate x
% y(13) IKs
% y(14) RyR (R)
% y(15) RyR (O)
% y(16) RyR (I)
% y(17) Na Buffer cleft
% y(18) Na Buffer sl
% y(19) Ca Buffer myofilament
% y(20) Ca Buffer TnCH
% y(21) Mg Buffer TnCH
% y(22) Ca Buffer CaM
% y(23) Ca Buffer Myosin
% y(24) Mg Buffer Myosin
% y(25) Ca Buffer SRB
% y(26) Ca Buffer - low cleft
% y(27) Ca Buffer - low sl
% y(28) Ca Buffer - high cleft
% y(29) Ca Buffer - high sl
% y(30) Ca-Csqn
% y(31) [Ca]sr
% y(32) [Na]cleft
% y(33) [Na]sl
% y(34) [Na]i
% y(35) [K]i (NOT USED)
% y(36) [Ca]cleft
% y(37) [Ca]sl
% y(38) [Ca]i
% y(39) Vm
% y(40) Itos - gate r
% y(41) INaL - gate m
% y(42) INaL - gate h

ydot = zeros(size(y));

%% Input
p_HR = p(1); % stimulation frequency (Hz)

%% Model Parameters
% Constants
R = 8314;       % [J/kmol*K]  
Frdy = 96485;   % [C/mol]  
Temp = 310;     % [K]
FoRT = Frdy/R/Temp;
Cmem = 1.3810e-10;   % [F] membrane capacitance
Qpow = (Temp-310)/10;

% Cell geometry
cellLength = 100;     % cell length [um]
cellRadius = 10.25;   % cell radius [um]
junctionLength = 160e-3;  % junc length [um]
junctionRadius = 15e-3;   % junc radius [um]
distSLcyto = 0.45;    % dist. SL to cytosol [um]
distJuncSL = 0.5;  % dist. junc to SL [um]
DcaJuncSL = 1.64e-6;  % Dca junc to SL [cm^2/sec]
DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
DnaJuncSL = 1.09e-5;  % Dna junc to SL [cm^2/sec]
DnaSLcyto = 1.79e-5;  % Dna SL to cyto [cm^2/sec] 
Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; Vsr = 0.035*Vcell; Vsl = 0.02*Vcell; Vjunc = 0.0539*.01*Vcell; 
SAjunc = 20150*pi*2*junctionLength*junctionRadius;  % [um^2]
SAsl = pi*2*cellRadius*cellLength;          % [um^2]
%J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 2.3056e-11
%J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 2.3056e-11
% tau's from c-code, not used here
J_ca_juncsl =1/1.2134e12; % [L/msec] = 8.2413e-13
J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] = 5.4621e-11

% Fractional currents in compartments
Fjunc = 0.11; Fsl = 1-Fjunc;
Fjunc_CaL = 0.9; Fsl_CaL = 1-Fjunc_CaL;

% Fixed ion concentrations     
Cli = 15;   % Intracellular Cl  [mM]
Clo = 150;  % Extracellular Cl  [mM]
Ko = 5.4;   % Extracellular K   [mM]
Nao = 140;  % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]
Mgi = 0.5;  % Intracellular Mg  [mM]

% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/y(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/y(33));       % [mV]
ek = (1/FoRT)*log(Ko/y(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/y(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/y(37));     % [mV]
ecl = (1/FoRT)*log(Cli/Clo);            % [mV]

%% Na transport parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag to choose Na channel
%flag=1 -> Markov
% flag=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNa =23;        % [mS/uF]
GNa = 16;
GNaL = 0.0045;
GNaB = 0.297e-3;    % [mS/uF] 
IbarNaK = 1.90719;     % [uA/uF]
KmNaip = 11;         % [mM]
KmKo = 1.5;         % [mM]
Q10NaK = 1.63;  
Q10KmNai = 1.39;

%% K current parameters
pNaK = 0.01833;      
GtoSlow = 0.06*1;     % [mS/uF] %0.09 CaMKII
GtoFast = 0.02*1;     % [mS/uF] 
gkp = 0.001;

%% Cl current parameters
GClCa = 0.109625;   % [mS/uF]
GClB = 9e-3;        % [mS/uF]
KdClCa = 100e-3;    % [mM]

%% Ca transport parameters
pCa = 1*5.4e-4;       % [cm/sec]
pNa = 1*1.5e-8;       % [cm/sec] 500* -> 7.5000e-006
%pNa = 3.6/4200*85/9*pCa; % 0.0081* -> 4.3714e-006
pK = 1*2.7e-7;        % [cm/sec]
Q10CaL = 1.8;       

IbarNCX = 1.0*9.0;      % [uA/uF]
KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]  
nu = 0.35;          % [none]
Kdact = 0.256e-3;   % [mM] 
Q10NCX = 1.57;      % [none]
IbarSLCaP = 0.0673; % IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
KmPCa = 0.5e-3;     % [mM] 
GCaB = 2.513e-4;    % [uA/uF] 
Q10SLCaP = 2.35;    % [none]

Q10SRCaP = 2.6;          % [none]
Vmax_SRCaP = 5.3114e-3;  % [mM/msec] (286 umol/L cytosol/sec)
Kmf = 0.246e-3;          % [mM] default
%Kmf = 0.175e-3;          % [mM]
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = 25;                 % [1/ms]      
koCa = 10;               % [mM^-2 1/ms]   %default 10   modified 20
kom = 0.06;              % [1/ms]     
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
ec50SR = 0.45;           % [mM]

%% Buffering parameters
Bmax_Naj = 7.561;       % [mM] % Bmax_Naj = 3.7; (c-code difference?)  % Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = 19.6e-3;    % [1/ms] 
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity 
koff_tnchca = 0.032e-3; % [1/ms] 
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms] 
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3;       % [mM] **? about setting to 0 in c-code**   % CaM buffering
koff_cam = 238e-3;      % [1/ms] 
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        % [mM]    % SL buffering
Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    % [mM]    %Fei *0.1!!! junction reduction factor
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       % [mM] 
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM] %Fei *0.1!!! junction reduction factor
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 140e-3*Vmyo/Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
koff_csqn = 65;         % [1/ms] 
kon_csqn = 100;         % [1/mM/ms] 

%% I_Na: Fast Na Current
am = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13)));
bm = 0.08*exp(-y(39)/11);
if y(39) >= -40,
    ah = 0; aj = 0;
    bh = 1/(0.13*(1+exp(-(y(39)+10.66)/11.1)));
    bj = 0.3*exp(-2.535e-7*y(39))/(1+exp(-0.1*(y(39)+32)));
else
    ah = 0.135*exp((80+y(39))/-6.8);
    bh = 3.56*exp(0.079*y(39))+3.1e5*exp(0.35*y(39));
    aj = (-1.2714e5*exp(0.2444*y(39))-3.474e-5*exp(-0.04391*y(39)))*(y(39)+37.78)/(1+exp(0.311*(y(39)+79.23)));
    bj = 0.1212*exp(-0.01052*y(39))/(1+exp(-0.1378*(y(39)+40.14)));
end
ydot(1) = am*(1-y(1))-bm*y(1);
ydot(2) = ah*(1-y(2))-bh*y(2);
ydot(3) = aj*(1-y(3))-bj*y(3);

%% I_NaL: Late Na Current
aml = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13)));
bml = 0.08*exp(-y(39)/11);
hlinf = 1/(1+exp((y(39)+91)/6.1));
tauhl = 600;
ydot(41) = aml*(1-y(41))-bml*y(41);
ydot(42) = (hlinf-y(42))/tauhl;

I_NaL_junc = Fjunc*GNaL*y(41)^3*y(42)*(y(39)-ena_junc);
I_NaL_sl = Fsl*GNaL*y(41)^3*y(42)*(y(39)-ena_sl);

%% INa MARKOV MODEL
% flag to switch between bGal and CaMKII
% bGal=1;
% 
% P1a1=3.802;
% P2a1=0.1027;
% % a2 
% P1a2=9.178;
% P2a2=25;
% % a3 
% P1a3= 3.7933*10^-7;
% P2a3=7.7;
% % b3 
% if bGal==1
% P1b3=.5*0.0084;
% else
% P1b3=.8*0.0084;
% end
% P2b3=.1*0.00002;
% % b1 
% P1b1=.1917;
% P2b1=20.3;
% % b12 
% P1b12=.2;
% P2b12=5;
% % b13 
% P1b13=.22;
% P2b13=10;
% % a4 b4 a5 b5 
% P1a4=100;
% if bGal==1
% P1a5=0.3543e-3;
% else
% P1a5=1.8*0.3543e-3;
% end
% P2a5=23.2696;
% if bGal==1
% P1b4=8.8061e-007;
% else
% P1b4=0.8*8.8061e-007;
% end
% P2b4=11.3944;
% P1b5=0.2868*10^-3;
% P2b5=35.9898;
% 
% 
% %Transition Rates
% vel=3;
% %C3 -> C2	       
% a11=vel*(P1a1/(P2a1*exp(-(y(39))/17)+0.20*exp(-(y(39))/150)));
% %C2 ->C1				
% a12=vel*(P1a1/(P2a1*exp(-(y(39))/15)+0.23*exp(-(y(39))/150)));
% %C1 ->O				
% a13=vel*(P1a1/(P2a1*exp(-(y(39))/12)+0.25*exp(-(y(39))/150)));
% 
% %C2 ->C3				
% b11=vel*P1b1*exp(-(y(39))/P2b1);
% %C1 -> C2		   	
% b12=vel*P1b12*exp(-(y(39)-P2b12)/(P2b1));
% %O -> C1	 			
% b13=vel*(P1b13*exp(-(y(39)-P2b13)/(P2b1))); 
% 
% %IC3 -> IC2	 		
% a111=a11;
% %IC2 -> IF	 	
% a112=a12;
% 
% %IC2 -> IC3   		
% b111=b11;
% %IF -> IC2
% b112=b12;
% 
% %IF -> C1    		
% a3=vel*P1a3*exp(-y(39)/P2a3);
% %IC2 ->C2   	   	
% a3=vel*P1a3*exp(-y(39)/P2a3);
% %IC3 ->C3   	   	
% a3=vel*P1a3*exp(-y(39)/P2a3);
% 
% %C1 -> IF				
% b3=vel*(P1b3+P2b3*y(39));
% %C2 -> IC2		
% b3=vel*(P1b3+P2b3*y(39));
% %C3 -> IC3			
% b3=vel*(P1b3+P2b3*y(39));
% 
% %O -> IF	
% a2=vel*(P1a2*exp(y(39)/P2a2));
% %IF ->O				
% b2=(a13*a2*a3)/(b13*b3);
% 
% %IF -> IM1			
% 
% a4=a2/P1a4;
% 
% %IM1 -> IF			
% b4=vel*P1b4*exp(-y(39)/P2b4);
% % b4=a3;
% %IM1 -> IM2			
% % a5=a2/P1a5;  
% a5=vel*P1a5*exp(y(39)/P2a5);
% 
% %IM2 -> IM1			
% b5=vel*P1b5*exp(-y(39)/P2b5);
% % b5=vel*P1b5/(1+exp((y(39)+20-P2b5)/P3b5));
% 
% % U->L
% if bGal==1
% a6=vel*4.7*1e-7; %default
% %a6=vel*1.6*4.7*1e-7;
% else
% a6=vel*2.5*4.7*1e-7; %default
% %a6=vel*3.2*4.7*1e-7;
% end
% % L->U
% b6=vel*9.5e-4;
% %2*9.5e-4;
% 
% %IM2=1-y(41)-y(42)-y(43)-y(44)-y(45)-y(46)-y(47)-y(48);
% %UIC3o; 41 UIC2o; 42 UIFo; 43 UIM1o; 44 UC3o; 45 UC2o; 46
% %UC1o; 47 UOo; 48 UIM2o; 49 LC3o; 50 LC2o; 51 LC1o; 52 LOo 53
% ydot(41) = -(a11+a3)*y(41) + b11*y(42) + b3*y(45); 
% ydot(42) = a11*y(41) - (a3+a12+b11)*y(42) + b12*y(43) + b3*y(46);
% ydot(43) = a12*y(42) - (b12+a4+b2+a3)*y(43) + b4*y(44) + b3*y(47) + a2*y(48);
% ydot(44) = a4*y(43) - (b4+a5)*y(44) + b5*y(49);
% ydot(45) = a3*y(41) - (b3+a11+a6)*y(45) + b11*y(46) + b6*y(50);
% ydot(46) = a3*y(42) + a11*y(45) - (b11+b3+a12+a6)*y(46) + b12*y(47) + b6*y(51);
% ydot(47) = a3*y(43) + a12*y(46) - (b3+b12+a13+a6)*y(47) + b13*y(48) + b6*y(52) ;
% ydot(48) = b2*y(43) + a13*y(47) - (a2+b13+a6)*y(48)  + b6*y(53);
% ydot(49) = a5*y(44) - b5*y(49);
% ydot(50) = -(b6+a11)*y(50) + a6*y(45) + b11*y(51);
% ydot(51) = -(b6+a12+b11)*y(51) + a6*y(46) + a11*y(50) + b12*y(52);
% ydot(52) = -(b6+a13+b12)*y(52) + a6*y(47) + a12*y(51) + b13*y(53);
% ydot(53) = -(b6+b13)*y(53) + a6*y(48) + a13*y(52);

% I_Na_junc2 = Fjunc*GNa*(y(48)+y(53))*(y(39)-ena_junc);
% I_Na_sl2 = Fsl*GNa*(y(48)+y(53))*(y(39)-ena_sl);

I_Na_junc = Fjunc*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_junc)+I_NaL_junc ;
I_Na_sl = Fsl*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_sl)+I_NaL_sl ;

% I_Na_junc= I_Na_junc1*(1-flag)+I_Na_junc2*flag;
% I_Na_sl= I_Na_sl1*(1-flag)+I_Na_sl2*flag;
I_Na = I_Na_junc+I_Na_sl;
I_NaL = I_NaL_junc+I_NaL_sl;
I_Natot = I_Na+I_NaL;

%% I_nabk: Na Background Current
I_nabk_junc = Fjunc*GNaB*(y(39)-ena_junc);
I_nabk_sl = Fsl*GNaB*(y(39)-ena_sl);
I_nabk = I_nabk_junc+I_nabk_sl;

%% I_nak: Na/K Pump Current
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0365*sigma*exp(-y(39)*FoRT));
I_nak_junc = Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/y(32))^4) /(Ko+KmKo);
I_nak_sl = Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/y(33))^4) /(Ko+KmKo);
I_nak = I_nak_junc+I_nak_sl;

%% I_kr: Rapidly Activating K Current
gkr = 1*0.03*sqrt(Ko/5.4);
xrss = 1/(1+exp(-(y(39)+50)/7.5));
tauxr = 1/(1.38e-3*(y(39)+7)/(1-exp(-0.123*(y(39)+7)))+6.1e-4*(y(39)+10)/(exp(0.145*(y(39)+10))-1));
ydot(12) = (xrss-y(12))/tauxr;
rkr = 1/(1+exp((y(39)+33)/22.4));
I_kr = gkr*y(12)*rkr*(y(39)-ek);

%% I_ks: Slowly Activating K Current
pcaks_junc = -log10(y(36))+3.0; 
pcaks_sl = -log10(y(37))+3.0;  
gks_junc = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
gks_sl = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6))); 
eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));	
xsss = 1/(1+exp(-(y(39)-1.5)/16.7));
tauxs = 1/(7.19e-5*(y(39)+30)/(1-exp(-0.148*(y(39)+30)))+1.31e-4*(y(39)+30)/(exp(0.0687*(y(39)+30))-1)); 
ydot(13) = (xsss-y(13))/tauxs;
I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);
I_ks = I_ks_junc+I_ks_sl;

%% I_kp: Plateau K current
kp_kp = 1/(1+exp(7.488-y(39)/5.98));
I_kp_junc = Fjunc*gkp*kp_kp*(y(39)-ek);
I_kp_sl = Fsl*gkp*kp_kp*(y(39)-ek);
I_kp = I_kp_junc+I_kp_sl;

%% I_to: Transient Outward K Current (slow and fast components)
xtoss = 1/(1+exp(-(y(39)+3.0)/15));
ytoss = 1/(1+exp((y(39)+33.5)/10));
rtoss = 1/(1+exp((y(39)+33.5)/10));
tauxtos = 9/(1+exp((y(39)+3.0)/15))+0.5;
tauytos = 3e3/(1+exp((y(39)+60.0)/10))+30;
%tauytos = 182/(1+exp((y(39)+33.5)/10))+1;
taurtos = 2.8e3/(1+exp((y(39)+60.0)/10))+220; %Fei changed here!! time-dependent gating variable
%taurtos = 8085/(1+exp((y(39)+33.5)/10))+313;
ydot(8) = (xtoss-y(8))/tauxtos;
ydot(9) = (ytoss-y(9))/tauytos;
ydot(40)= (rtoss-y(40))/taurtos; %Fei changed here!! time-dependent gating variable
I_tos = GtoSlow*y(8)*(y(9)+0.5*y(40))*(y(39)-ek); % [uA/uF]

tauxtof = 3.5*exp(-y(39)*y(39)/30/30)+1.5;
%tauxtof = 3.5*exp(-((y(39)+3)/30)^2)+1.5;
tauytof = 20.0/(1+exp((y(39)+33.5)/10))+20.0;
%tauytof = 20.0/(1+exp((y(39)+33.5)/10))+20.0;
ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;
I_tof = GtoFast*y(10)*y(11)*(y(39)-ek);

I_to = I_tos + I_tof;

%% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(y(39)-ek-59.215)));
bki =(0.49124*exp(0.08032*(y(39)+5.476-ek)) + exp(0.06175*(y(39)-ek-594.31))) /(1 + exp(-0.5143*(y(39)-ek+4.753)));
kiss = aki/(aki+bki);
I_ki = 1*0.9*sqrt(Ko/5.4)*kiss*(y(39)-ek);

%% I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y(36))*(y(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y(37))*(y(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;
I_Clbk = GClB*(y(39)-ecl);

%% I_Ca: L-type Calcium Current
dss = 1/(1+exp(-(y(39)+14.5)/6.0));
taud = dss*(1-exp(-(y(39)+14.5)/6.0))/(0.035*(y(39)+14.5));
fss = 1/(1+exp((y(39)+35.06)/3.6))+0.6/(1+exp((50-y(39))/20));
tauf = 1/(0.0197*exp( -(0.0337*(y(39)+14.5))^2 )+0.02);
ydot(4) = (dss-y(4))/taud;
ydot(5) = (fss-y(5))/tauf;
ydot(6) = 1.7*y(36)*(1-y(6))-11.9e-3*y(6); % fCa_junc   koff!!!!!!!!
ydot(7) = 1.7*y(37)*(1-y(7))-11.9e-3*y(7); % fCa_sl
fcaCaMSL = 0.1/(1+(0.01/y(37)));
fcaCaj = 0.1/(1+(0.01/y(36)));
fcaCaMSL =0;
fcaCaj = 0;
ibarca_j = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibarca_sl = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibark = pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
ibarna_j = pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
ibarna_sl = pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
I_Ca_junc = (Fjunc_CaL*ibarca_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45*1;
I_Ca_sl = (Fsl_CaL*ibarca_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.45*1;
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK = (ibark*y(4)*y(5)*(Fjunc_CaL*(fcaCaj+(1-y(6)))+Fsl_CaL*(fcaCaMSL+(1-y(7))))*Q10CaL^Qpow)*0.45*1;
I_CaNa_junc = (Fjunc_CaL*ibarna_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45*1;
I_CaNa_sl = (Fsl_CaL*ibarna_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*.45*1;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca+I_CaK+I_CaNa;

%% I_ncx: Na/Ca Exchanger flux
Ka_junc = 1/(1+(Kdact/y(36))^3);
Ka_sl = 1/(1+(Kdact/y(37))^3);
s1_junc = exp(nu*y(39)*FoRT)*y(32)^3*Cao;
s1_sl = exp(nu*y(39)*FoRT)*y(33)^3*Cao;
s2_junc = exp((nu-1)*y(39)*FoRT)*Nao^3*y(36);
s3_junc = KmCai*Nao^3*(1+(y(32)/KmNai)^3) + KmNao^3*y(36)*(1+y(36)/KmCai)+KmCao*y(32)^3+y(32)^3*Cao+Nao^3*y(36);
s2_sl = exp((nu-1)*y(39)*FoRT)*Nao^3*y(37);
s3_sl = KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37)*(1+y(37)/KmCai)+KmCao*y(33)^3+y(33)^3*Cao+Nao^3*y(37);
I_ncx_junc = Fjunc*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx_sl = Fsl*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx = I_ncx_junc+I_ncx_sl;

%% I_pca: Sarcolemmal Ca Pump Current
I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y(36)^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y(37)^1.6);
I_pca = I_pca_junc+I_pca_sl;

%% I_cabk: Ca Background Current
I_cabk_junc = Fjunc*GCaB*(y(39)-eca_junc);
I_cabk_sl = Fsl*GCaB*(y(39)-eca_sl);
I_cabk = I_cabk_junc+I_cabk_sl;

%% SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y(31))^2.5);
koSRCa = koCa/kCaSR;
kiSRCa = kiCa*kCaSR;
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-kom*y(15));   % R
ydot(15) = (koSRCa*y(36)^2*y(14)-kom*y(15))-(kiSRCa*y(36)*y(15)-kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-kim*y(16))-(kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = ks*y(15)*(y(31)-y(36));          % [mM/ms]

J_serca = Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    /(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP);
J_SRleak = 5.348e-6*(y(31)-y(36));           %   [mM/ms]

%% Sodium and Calcium Buffering
ydot(17) = kon_na*y(32)*(Bmax_Naj-y(17))-koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*y(33)*(Bmax_Nasl-y(18))-koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = kon_tncl*y(38)*(Bmax_TnClow-y(19))-koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(20) = kon_tnchca*y(38)*(Bmax_TnChigh-y(20)-y(21))-koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-y(20)-y(21))-koff_tnchmg*y(21);   % TnCHm     [mM/ms]
ydot(22) = kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*y(38)*(Bmax_myosin-y(23)-y(24))-koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-y(23)-y(24))-koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*y(38)*(Bmax_SR-y(25))-koff_sr*y(25);                    % SRB       [mM/ms]
%J_CaB_cytosol = sum(ydot(19:25)); % wrong formulation
J_CaB_cytosol = ydot(19)+ydot(20)+ydot(22)+ydot(23)+ydot(25);

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*y(36)*(Bmax_SLlowj-y(26))-koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*y(37)*(Bmax_SLlowsl-y(27))-koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*y(36)*(Bmax_SLhighj-y(28))-koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*y(37)*(Bmax_SLhighsl-y(29))-koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

%% Ion concentrations
% SR Ca Concentrations
ydot(30) = kon_csqn*y(31)*(Bmax_Csqn-y(30))-koff_csqn*y(30);       % Csqn      [mM/ms]
ydot(31) = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30);         % Ca_sr     [mM/ms] %Ratio 3 leak current

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
%I_Na_tot_sl2 = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl*0;   % [uA/uF]

ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y(32)-y(33))...
   +J_na_slmyo/Vsl*(y(34)-y(33))-ydot(18);
ydot(34) = J_na_slmyo/Vmyo*(y(33)-y(34));             % [mM/msec] 

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp;     % [uA/uF]
ydot(35) = 0; %-I_K_tot*Cmem/(Vmyo*Frdy);           % [mM/msec]

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y(37)-y(36))...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y(36)-y(37))...
    + J_ca_slmyo/Vsl*(y(38)-y(37))-J_CaB_sl;   % Ca_sl
ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol+J_ca_slmyo/Vmyo*(y(37)-y(38));

%% Simulation type
protocol = 'pace';

switch lower(protocol)
    case {'none',''},
        I_app = 0;

    case 'pace', % pace w/ current injection at rate 'rate'
		rate = (p_HR)*1e-3;
		if mod(t+0,1/rate) <= 5
            I_app = 9.5;
        else
            I_app = 0.0;
        end
        
    case 'vclamp',      
		V_hold = -90;
        V_test = 80;
		if (t > 100 && t < 2100)
		    V_clamp = V_test;
		else
		    V_clamp = V_hold;
		end
		R_clamp = 0.02;
		I_app = (V_clamp-y(39))/R_clamp;
end 

%% Membrane Potential
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk;                        % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
ydot(39) = -(I_tot-I_app);
vmax = ydot(39);
% ----- END EC COUPLING MODEL ---------------

% adjust output depending on the function call
if (nargin == 3)
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'ydot')
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'rates')
    output = r;
elseif (nargin == 4) && strcmp(runType,'currents')
    %currents = [I_Na I_nabk I_nak I_kr I_ks I_kp I_tos I_tof I_ki I_ClCa I_Clbk I_Catot I_ncx I_pca I_cabk J_serca*Vmyo/Vsr];
    currents = [I_Catot I_Natot vmax];
    output = currents;
end

%% Calculate timecourse for currents and other intermediates
function currents = calcCurrents(t,y,p)
% After running a simulation, feed the time vector and state variables into
% this function to compute ionic currents, etc.
% currents = [I_Catot I_Natot];
currents = [];
for i=1:size(t)
    if ceil(i/1000)==i/1000
        disp(['t = ',num2str(ceil(t(i)))]);
    end
    currents = [currents;f(t(i),y(i,:),p,'currents')];
end
% end calcCurrents