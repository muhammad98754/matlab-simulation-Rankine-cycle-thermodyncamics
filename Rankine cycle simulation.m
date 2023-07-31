% Program to simulate Rankine cycle in matlab for interview preparation at Pepper motion gmbh

%-----------------------------------------------------
%References

% Engineering Thermodynamics,PI Dhar : Concepts of Thermodynamics
% ENGINEERING THERMODYNAMICS, R. K. Rajput :VAPOUR POWER CYCLES
% Wikipedia : Rankine cycle


%-----------------------------------------------------
clear all
close all
clc

fprintf('%50s','Rankine Cycle Simulation')

%-----------------------------------------------------
%Rankine cycle satges


fprintf('\n\n stage 1 - 2 is Isentropic Expansion in the Turbine')
fprintf('\n stage 2 - 3 Constant Pressure Heat Rejection by the Condenser')
fprintf('\n stage 3 - 4 is Isentropic Compression in the Pump')
fprintf('\n stage 4 - 1 is Constant Pressure Heat Addition by the Boiler')

%-----------------------------------------------------
%user input

p1 = 30  % Pressure at the Turbine Inlet (in bar)
t1 = 400 % Temperature at the Turbine Inlet (in Degree Celcius)
p2 = 0.05 % Pressure at the Condenser Inlet (in bar)

%-----------------------------------------------------
% from XSteam library calculate enthalpy and entropy values at different stages

% at turbine inlet
s1 = XSteam('s_pT',p1,t1)
h1 = XSteam('h_pT',p1,t1)

% at Condenser Inlet
s2 = s1;   % stage 1-2 Isentropic Expansion process in the Turbine
sf2=XSteam('sL_p',p2);  %entropy of liquid phase at given pressure p2
sg2=XSteam('sV_p',p2);   %entropy of vapor phase at given pressure p2
hf2=XSteam('hL_p',p2);   %enthalpy of liquid phase at given pressure p2
hg2=XSteam('hV_p',p2);   %enthalpy of vapor phase at given pressure p2

%`Determining Dryness fraction`  [The dryness fraction means 
% the amount of vapor present in a two-phase mixture.
% It is a value between 0 and 1, where 0 indicates a completely liquid state
% and 1 indicates a completely vapor state.]

if(s2 <= sg2)
    x2 = ((s2 - sf2) / (sg2 - sf2));
else
    x2 = 1;    
end
h2 = hf2 + (x2 * (hg2 - hf2));

%Pump Inlet
p3 = p2;  % Constant Pressure Heat Rejection by the Condenser
t3 = XSteam('Tsat_p',p3);
s3 = XSteam('sL_p',p3);
h3 = XSteam('hL_p',p3);
t2 = t3;  % From the T-s diagram, it is clear that t3 = t2.

%Boiler Inlet 
p4 = p1; 
s4 = s3; 
h4 = XSteam('h_ps',p4,s4);
%amount of heat required to change the temperature of a substance
% at constant pressure (cp) or constant volume (cv)
cp = XSteam('Cp_ps',p4,s4);   %Isobaric specific heat
cv = XSteam('Cv_ps',p4,s4);   %Isochoric specific heat
%Determining Gamma
% Gamma means how the temperature and pressure of a substance 
% are related when no heat is exchanged with the surroundings.
gamma = ( cp / cv );
t4 = (t3 * (( p3 / p4 )^(( 1 - gamma) / gamma)));

%Work Done by the Turbine
Wt = h1 - h2;
%Work Done by the Pump
Wp = h4 - h3;
%Net Work Done
Wnet = Wt - Wp;

%Heat Added in the Boiler
qin = h1 - h4;
%Heat Rejected in the Condensor
qout = h2 - h3;

%Thermal Efficiency
n = (Wnet / qin) * 100;

%Specific Steam Consumption
SSC = (3600/Wnet);  %how much steam is needed to generate
% a certain amount of work 



%-----------------------------------------------------
%results

fprintf('\n %45s ','Result');
fprintf('\n %s %50s','At State Point 1','At State Point 4');
fprintf('\n%s %.2f%s %40s %.2f%s',' P1 is :',p1,' Bar',' P4 is :',p4,' Bar');
fprintf('\n%s %.2f%s %32s %.2f%s',' T1 is :',t1,' Deg Celcius',' T4 is :',t4,' Deg Celcius');
fprintf('\n%s %.2f%s %38s %.2f%s',' H1 is :',h1,' kJ/kg',' H4 is :',h4,' kJ/kg');
fprintf('\n%s %.2f%s %37s %.2f%s',' S1 is :',s1,' kJ/kgK',' S4 is :',s4,' kJ/kgK');


fprintf('\n%s %42s %.2f%s','At State Point 2','Wt is :',Wt,' kJ/kg');
fprintf('\n%s %.2f%s %40s %.2f%s',' P2 is :',p2,' Bar','Wp is :',Wp,' kJ/kg');
fprintf('\n%s %.2f%s %33s %.2f%s',' T2 is :',t2,' Deg Celcius','Wnet is :',Wnet,' kJ/kg');
fprintf('\n%s %.2f%s %41s %.2f%s',' H2 is :',h2,' kJ/kg','Ntherm is :',n,' Percent');
fprintf('\n%s %.2f%s %41s %.2f%s',' S2 is :',s2,' kJ/kgK','S.S.C is :',SSC,' kg/kWh');
fprintf('\n%8s %.2f%s','X2 is :',x2);

fprintf('\n%s','At State Point 4');
fprintf('\n%s %.2f%s',' P1 is :',p4,' Bar');
fprintf('\n%s %.2f%s',' T1 is :',t4,' Deg Celcius');
fprintf('\n%s %.2f%s',' H1 is :',h4,' kJ/kg');
fprintf('\n%s %.2f%s',' S1 is :',s4,' kJ/kgK');


%-----------------------------------------------------
% T-s diagram plot

t = linspace(1,t1,2000);  %temperature is converted to an array of temperatures
%find liquid entropy (sl) and vapor entropy (sv) at every point of the temperature
for i = 1: 2000
    sl(i) = XSteam('sL_T',t(i));
end
for i = 1: 2000
    sv(i) = XSteam('sV_T',t(i));
end 

figure(1)
plot(sl,t,'-.r','LineWidth',2)
hold on
plot(sv,t,'-.r','LineWidth',2)
title(' T-s Diagram of an Ideal Rankine Cycle')
xlabel('Entropy(KJ/KgK)')
ylabel('Temperature(k)')

%1-2 State
plot([s1 s2],[t1 t2],'k','LineWidth',1)
%2-3 State
plot([s2 s3],[t2 t3],'k','LineWidth',1)
%3-4 State
plot([s3 s4],[t3 t4],'k','LineWidth',1)
%4-1 State
tp = XSteam('TSat_p',p4); %calculates the saturation temperature 
sp = XSteam('sL_T',tp); %entropy of the saturated liquid phase
tq = tp;
sq = XSteam('sV_T',tq); %entropy of the saturated vapor phase
plot([s4,sp],[t4,tp],'k','LineWidth',1)
plot([sp sq],[tp tq],'k','LineWidth',1)
plot([sq,s1],[tq,t1],'k','LineWidth',1)
%name the points
text(s1,t1,' 1');
text(s2,t2,' 2');
text(s3,t3,' 3');
text(s4,t4,' 4');



%-----------------------------------------------------
%h-s diagram  plot
%for this determine liquid enthalpy (hl) and vapor enthalpy (hv)

for i = 1: 2000
    hl(i) = XSteam('hL_T',t(i)); %enthalpy of the saturated liquid phase
end
for i = 1: 2000
    hv(i) = XSteam('hV_T',t(i)); %enthalpy of the saturated vapor phase 
end
figure(2)

plot(sl,hl,'-.r','LineWidth',2)
hold on
plot(sv,hv,'-.r','LineWidth',2)
title(' H-S Diagram of an Ideal Rankine Cycle')
xlabel('Entropy(KJ/KgK)')
ylabel('Enthalpy(KJ/Kg)')

% plot state 1 to state 2, state 2 to state 3, state 3 to state 4, and state 4 to state 1.
%1-2 State
plot([s1 s2],[h1 h2],'k','LineWidth',1)
%2-3 State
plot([s2 s3],[h2 h3],'k','LineWidth',1)
%3-4 State
plot([s3 s4],[h3 h4],'k','LineWidth',1)
%4-1 State
plot([s4 s1],[h4 h1],'k','LineWidth',1)

%name the points
text(s1,h1,' 1');
text(s2,h2,' 2');
text(s3,h3,' 3');
text(s4,h4,' 4');