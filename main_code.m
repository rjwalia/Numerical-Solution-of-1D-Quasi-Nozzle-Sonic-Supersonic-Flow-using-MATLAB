clc
clear all 
close all

%Rajat Walia
%SOLVING QUASI 1D SUPERSONIC NOZZLE FLOW
%All flow field variables are calculated and plotted in their non-dimensional form

%Inputs

n = 31;                %mesh points
x = linspace(0,3,n);   %1d mesh
dx = x(2) - x(1);      %grid spacing

gamma = 1.4;           %gas constant

nt  = 800;             %number of time step
CFL = 0.5;             %CFL No.

% NON - CONSERVATIVE FORM

[v,pressure,T,Mach_no,mass_flow,rho,rho_throat,P_throat,v_throat,T_throat,Mach_no_throat, mass_flow_throat] = non_conservative_form(n,x,dx,gamma,CFL,nt);

%Potting time wise variation of non-dimensional primitive variables at the nozzle throat

hold on
figure(4)

%density
subplot(4,1,1)
plot(linspace(1,nt,nt),rho_throat,'color','m', 'LineWidth', 1.5)
ylabel('Density')
axis([0 1000 0 1.4]);
grid minor;
title("Non Dimensional Flow Variable variation with Time at Nozzle's Throat [Non Conservative Form]")

%temperature
subplot(4,1,2)
plot(linspace(1,nt,nt),T_throat,'color','c', 'LineWidth', 1.5)
ylabel('Temperature')
axis([0 1000 0.5 1]);
grid minor;

%pressure
subplot(4,1,3)
plot(linspace(1,nt,nt),P_throat,'color','g', 'LineWidth', 1.5')
ylabel('Pressure')
axis([0 1000 0.3 1.2]);
grid minor;

%Mach number
subplot(4,1,4)
plot(linspace(1,nt,nt),Mach_no_throat,'color','y', 'LineWidth', 1.5)
xlabel('Time Step')
ylabel('Mach No.')
axis([0 1000 0.5 1.5]);
grid minor;

%Potting steady-state field of non-dimensional primitive variables across nozzle

hold on
figure(5)

%density
subplot(4,1,1)
plot(x,rho,'color','m', 'LineWidth', 1.5)
ylabel('Density')
axis([0 3 0 1]);
grid minor;
title("Steady-State Flow Field vs Nozzle X-Direction [Non Conservative Form]")

%temperature
subplot(4,1,2)
plot(x,T,'color','c', 'LineWidth', 1.5)
ylabel('Temperature')
axis([0 3 0 1]);
grid minor;

%pressure
subplot(4,1,3)
plot(x,pressure,'color','g', 'LineWidth', 1.5')
ylabel('Pressure')
axis([0 3 0 1]);
grid minor;

%Mach number
subplot(4,1,4)
plot(x,Mach_no,'color','y', 'LineWidth', 1.5)
xlabel('Nozzle X-Direction')
ylabel('Mach No.')
axis([0 3 0 4]);
grid minor;

mass_flow_rate_NC  = mass_flow;
T_NC    = T_throat;
rho_NC  = rho_throat;
P_NC    = P_throat;
Mach_NC = Mach_no_throat;

%CONSERVATIVE - FORM

[v,pressure,T,Mach_no,mass_flow,rho,rho_throat,P_throat,v_throat,T_throat,Mach_no_throat, mass_flow_throat] = conservative_form(n,x,dx,gamma,CFL,nt);

%Potting time wise variation of non-dimensional primitive variables at the nozzle throat

hold on
figure(7)

%density
subplot(4,1,1)
plot(linspace(1,nt,nt),rho_throat,'color','m', 'LineWidth', 1.5)
ylabel('Density')
axis([0 1000 0 1.4]);
grid minor;
title("Non Dimensional Flow Variable variation with Time at Nozzle's Throat [Conservative Form]")

%temperature
subplot(4,1,2)
plot(linspace(1,nt,nt),T_throat,'color','c', 'LineWidth', 1.5)
ylabel('Temperature')
axis([0 1000 0.5 1]);
grid minor;

%pressure
subplot(4,1,3)
plot(linspace(1,nt,nt),P_throat,'color','g', 'LineWidth', 1.5')
ylabel('Pressure')
axis([0 1000 0.3 1.2]);
grid minor;

%Mach number
subplot(4,1,4)
plot(linspace(1,nt,nt),Mach_no_throat,'color','y', 'LineWidth', 1.5)
xlabel('Time Step')
ylabel('Mach No.')
axis([0 1000 0.5 1.5]);
grid minor;

%Potting steady-state field of non-dimensional primitive variables across the nozzle

hold on
figure(8)

%density
subplot(4,1,1)
plot(x,rho,'color','m', 'LineWidth', 1.5)
ylabel('Density')
axis([0 3 0 1]);
grid minor;
title("Steady-State Flow Field vs Nozzle X-Direction [Conservative Form]")

%temperature
subplot(4,1,2)
plot(x,T,'color','c', 'LineWidth', 1.5)
ylabel('Temperature')
axis([0 3 0 1]);
grid minor;

%pressure
subplot(4,1,3)
plot(x,pressure,'color','g', 'LineWidth', 1.5')
ylabel('Pressure')
axis([0 3 0 1]);
grid minor;

%Mach number
subplot(4,1,4)
plot(x,Mach_no,'color','y', 'LineWidth', 1.5)
xlabel('Nozzle X-Direction')
ylabel('Mach No.')
axis([0 3 0 4]);
grid minor;

mass_flow_rate_C  = mass_flow;
T_C    = T_throat;
rho_C  = rho_throat;
P_C    = P_throat;
Mach_C = Mach_no_throat;

%Comparison of normalized steady-state mass-flow rate across the nozzle 
%between Conservative & Non-Conservative Form

%analytical non dimensional mass flow rate throught the nozzle =  0.590
mass_flow_rate_analytical = 0.590*ones(1,n); 

figure(9)
hold on
plot(x, mass_flow_rate_NC, 'r', 'LineWidth', 1.5); 
hold on
plot(x, mass_flow_rate_C, 'k', 'LineWidth', 1.5); 
hold on
plot(x, mass_flow_rate_analytical, 'g', 'LineWidth', 1.5, 'LineStyle','--'); 
grid minor;

legend("Non Conservative Form", "Conservative Form", "Analytical");
xlabel("Nozzle X-Direction");
ylabel("Non Dimensional Mass Flow Rate");
title("Comparison of steady Normalized Mass Flow Rate for Conservative & Non-Conservative Form")