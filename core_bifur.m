function [x,v,s,h,f] = core_bifur

curdir = pwd;
init;
cd(curdir);
opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',500000);
opt=contset(opt,'MinStepsize',1);
opt=contset(opt,'MaxStepsize',100);
opt=contset(opt,'Eigenvalues',1);

% Degradation rate:
ks = 0.125;  ku200 = 0.05;   kmz = 0.5;   kz = 0.1;
kmsl=0.5; ksl=0.1155; kk=0.1732;
% Transcription rate:
gs = 18000; gu200 = 2100;   gmz = 11;   gz = 100; 
gmsl=90; gsl=50000; gk=50000;
% Hills function threshold :
I0s=100000; z0u200 = 220000;   z0mz = 27500;   s0u200 = 180000;   s0mz = 180000; u2000 = 10000; sl0u200=220000; 
sl0msl=150000; sl0s=225000; s0msl=180000; s0s=300000; k0s=275000; k0msl=300000; s0k=180000; sl0k=225000; k0k=250000;
% Cooperativity:
nsmz = 2;  nIs = 2;  nzu200 = 3;   nsu200 = 2;   nzmz = 2;   nu200 = 6;  
nslu200=1; nslmsl=4; nsls=3; nsmsl=1; nss=5; nks=2; nsk=2; nslk=4; nkmsl=2; nkk=3;
% fold change
lamdazu200 =0.1;   lamdasu200 = 0.1;  lamdazmz = 7.5;   lamdasmz = 10; lamdaIs=3;
lamdaslu200=0.4; lamdaslmsl=2; lamdasls=0.5; lamdasmsl=0.5; lamdass=0.4;
lamdakk=2; lamdask=0.25; lamdaslk=0.5; lamdakmsl=0.25; lamdaks=0.5;

% external signal
I=0;



ap = 1; %describes the index of parameter for which the bifurcation is drawn using the init_EP_EP function. Currently, ap=1, thus bifurcation parameter is s (SNAIL levels)
handles = feval(@core);
tspan = 0:100:50000;

% initial condition
x_start = [33554.833280 56.5 0 0 0 0 0];

%calculating steady state for given initial condition 
[t,x_time] = ode15s(@(t,kmrgd)handles{2}(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k),tspan,x_start);
x_init = x_time(end,:)';

%drawing bifurcation using a continuation method
[x0,v0] = init_EP_EP(@core,x_init,[I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k],ap);
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);