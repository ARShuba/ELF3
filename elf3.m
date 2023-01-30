function out = elf3
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
%-----------------------------------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k,ge,ke,lamdaee,e0e,nee,lamdasle,sl0e,nsle,lamdase,s0e,nse,lamdaemz,e0mz,nemz)

%Ym function components
Mu0=1/(1+kmrgd(1)/u2000)^nu200;
Mu1=(kmrgd(1)/u2000)/(1+kmrgd(1)/u2000)^nu200;
Mu2=(kmrgd(1)/u2000)^2/(1+kmrgd(1)/u2000)^nu200;
Mu3=(kmrgd(1)/u2000)^3/(1+kmrgd(1)/u2000)^nu200;
Mu4=(kmrgd(1)/u2000)^4/(1+kmrgd(1)/u2000)^nu200;
Mu5=(kmrgd(1)/u2000)^5/(1+kmrgd(1)/u2000)^nu200;
Mu6=(kmrgd(1)/u2000)^6/(1+kmrgd(1)/u2000)^nu200;

%Hills functions
Hillszu200=(1+lamdazu200*(kmrgd(3)/z0u200)^nzu200)/(1+(kmrgd(3)/z0u200)^nzu200);
Hillssu200=(1+lamdasu200*(kmrgd(4)/s0u200)^nsu200)/(1+(kmrgd(4)/s0u200)^nsu200);
Hillszmz=(1+lamdazmz*(kmrgd(3)/z0mz)^nzmz)/(1+(kmrgd(3)/z0mz)^nzmz);
Hillssmz=(1+lamdasmz*(kmrgd(4)/s0mz)^nsmz)/(1+(kmrgd(4)/s0mz)^nsmz);
HillsIs=(1+lamdaIs*(I/I0s)^nIs)/(1+(I/I0s)^nIs);
Hillsss=(1+lamdass*(kmrgd(4)/s0s)^nss)/(1+(kmrgd(4)/s0s)^nss);
Hillssls=(1+lamdasls*(kmrgd(6)/sl0s)^nsls)/(1+(kmrgd(6)/sl0s)^nsls);
Hillssmsl=(1+lamdasmsl*(kmrgd(4)/s0msl)^nsmsl)/(1+(kmrgd(4)/s0msl)^nsmsl);
Hillsslu200=(1+lamdaslu200*(kmrgd(6)/sl0u200)^nslu200)/(1+(kmrgd(6)/sl0u200)^nslu200);
Hillsslk=(1+lamdaslk*(kmrgd(6)/sl0k)^nslk)/(1+(kmrgd(6)/sl0k)^nslk);
Hillssk=(1+lamdask*(kmrgd(4)/s0k)^nsk)/(1+(kmrgd(4)/s0k)^nsk);
Hillsks=(1+lamdaks*(kmrgd(7)/k0s)^nks)/(1+(kmrgd(7)/k0s)^nks);
Hillskmsl=(1+lamdakmsl*(kmrgd(7)/k0msl)^nkmsl)/(1+(kmrgd(7)/k0msl)^nkmsl);
Hillskk=(1+lamdakk*(kmrgd(7)/k0k)^nkk)/(1+(kmrgd(7)/k0k)^nkk);
Hillsemz = (1+lamdaemz*(kmrgd(8)/e0mz)^nemz)/(1+(kmrgd(8)/e0mz)^nemz);
Hillsse = (1+lamdase*(kmrgd(4)/s0e)^nse)/(1+(kmrgd(4)/s0e)^nse);
Hillssle = (1+lamdasle*(kmrgd(6)/sl0e)^nsle)/(1+(kmrgd(6)/sl0e)^nsle);
Hillsee = (1+lamdaee*(kmrgd(8)/e0e)^nee)/(1+(kmrgd(8)/e0e)^nee);

dydt=[gu200*Hillszu200*Hillssu200*Hillsslu200-kmrgd(2)*(0.005*6*Mu1+2*0.05*15*Mu2+3*0.5*20*Mu3+4*0.5*15*Mu4+5*0.5*6*Mu5+6*0.5*Mu6)-kmrgd(5)*(0.005*6*Mu1+2*0.05*15*Mu2+3*0.5*20*Mu3+4*0.5*15*Mu4+5*0.5*6*Mu5+6*0.5*Mu6)-ku200*kmrgd(1);
gmz*Hillszmz*Hillssmz*Hillsemz-kmrgd(2)*(0.04*6*Mu1+0.2*15*Mu2+20*Mu3+15*Mu4+6*Mu5+Mu6)-kmz*kmrgd(2);
gz*kmrgd(2)*(Mu0+0.6*6*Mu1+0.3*15*Mu2+0.1*20*Mu3+0.05*15*Mu4+0.05*6*Mu5+0.05*Mu6)-kz*kmrgd(3)
gs*HillsIs*Hillsss*Hillssls*Hillsks-ks*kmrgd(4);
gmsl*Hillssmsl*Hillskmsl-kmrgd(5)*(0.04*6*Mu1+0.2*15*Mu2+20*Mu3+15*Mu4+6*Mu5+Mu6)-kmsl*kmrgd(5);
gsl*kmrgd(5)*(Mu0+0.6*6*Mu1+0.3*15*Mu2+0.1*20*Mu3+0.05*15*Mu4+0.05*6*Mu5+0.05*Mu6)-ksl*kmrgd(6);
gk*Hillsslk*Hillssk*Hillskk-kk*kmrgd(7);
ge*Hillsse*Hillssle*Hillsee-ke*kmrgd(8)
] ;
    
    % --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(elf3);
y0=[0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k,ge,ke,lamdaee,e0e,nee,lamdasle,sl0e,nsle,lamdase,s0e,nse,lamdaemz,e0mz,nemz)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k,ge,ke,lamdaee,e0e,nee,lamdasle,sl0e,nsle,lamdase,s0e,nse,lamdaemz,e0mz,nemz)
%--------------------------------------------------------------------------
function hess = hessians(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k,ge,ke,lamdaee,e0e,nee,lamdasle,sl0e,nsle,lamdase,s0e,nse,lamdaemz,e0mz,nemz)
%--------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k,ge,ke,lamdaee,e0e,nee,lamdasle,sl0e,nsle,lamdase,s0e,nse,lamdaemz,e0mz,nemz)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k,ge,ke,lamdaee,e0e,nee,lamdasle,sl0e,nsle,lamdase,s0e,nse,lamdaemz,e0mz,nemz)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k,ge,ke,lamdaee,e0e,nee,lamdasle,sl0e,nsle,lamdase,s0e,nse,lamdaemz,e0mz,nemz)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdaslmsl,nslmsl,sl0msl,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl,k0msl,lamdakmsl,nkmsl,k0s,lamdaks,nks,s0k,lamdask,nsk,sl0k,lamdaslk,nslk,gk,kk,lamdakk,nkk,k0k,ge,ke,lamdaee,e0e,nee,lamdasle,sl0e,nsle,lamdase,s0e,nse,lamdaemz,e0mz,nemz)