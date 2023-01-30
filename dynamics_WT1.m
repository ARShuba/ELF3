function dydt = dynamics_WT1(t,kmrgd)

% Degradation rate:
ks = 0.125;  ku200 = 0.05;   kmz = 0.5;   kz = 0.1;
kmsl=0.5; ksl=0.1155; kk=0.1732; ke=0.125; kw=0.2722;
% Transcription rate:
gs = 18000; gu200 = 2100;   gmz = 11;   gz = 100; 
gmsl=90; gsl=50000; gk=50000; ge=50000; gw=10000;
% Hills function threshold :
I0s=100000; z0u200 = 220000;   z0mz = 27500;   s0u200 = 180000;   s0mz = 180000; u2000 = 10000; sl0u200=220000; 
sl0msl=150000; sl0s=225000; s0msl=180000; s0s=300000; k0s=275000; k0msl=300000; s0k=180000; sl0k=250000; k0k=275000;
e0e=200000; sl0e=220000; s0e=250000; e0mz=180000; w0s=350000; w0msl=220000; w0w=280000;
% Cooperativity:
nsmz = 2;  nIs = 2;  nzu200 = 3;   nsu200 = 2;   nzmz = 2;   nu200 = 6;  
nslu200=1; nslmsl=4; nsls=3; nsmsl=1; nss=5; nks=2; nsk=2; nslk=4; nkmsl=2; nkk=3;
nee=3; nsle=4; nse=4; nemz=2; nww=8; nwmsl=1; nws=1;
% fold change
lamdazu200 =0.1;   lamdasu200 = 0.1;  lamdazmz = 7.5;   lamdasmz = 10; lamdaIs=3;
lamdaslu200=0.4; lamdaslmsl=2; lamdasls=0.5; lamdasmsl=0.5; lamdass=0.4;
lamdakk=2; lamdask=0.25; lamdaslk=0.5; lamdakmsl=0.25; lamdaks=0.5; 
lamdaee=4; lamdasle=0.3; lamdase=0.25; lamdaemz=0.675; lamdaws=6; lamdawmsl=0.5; lamdaww=0.4;

% external signal
if t<=200
    I=0;
elseif t>1300
    I=0;
else
    I=160000;
end

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
Hillsws =(1+lamdaws*(kmrgd(9)/w0s)^nws)/(1+(kmrgd(9)/w0s)^nws);
Hillswmsl =(1+lamdawmsl*(kmrgd(9)/w0msl)^nwmsl)/(1+(kmrgd(9)/w0s)^nwmsl);
Hillsww =(1+lamdaww*(kmrgd(9)/w0w)^nww)/(1+(kmrgd(9)/w0w)^nww);

dydt=[gu200*Hillszu200*Hillssu200*Hillsslu200-kmrgd(2)*(0.005*6*Mu1+2*0.05*15*Mu2+3*0.5*20*Mu3+4*0.5*15*Mu4+5*0.5*6*Mu5+6*0.5*Mu6)-kmrgd(5)*(0.005*6*Mu1+2*0.05*15*Mu2+3*0.5*20*Mu3+4*0.5*15*Mu4+5*0.5*6*Mu5+6*0.5*Mu6)-ku200*kmrgd(1);
gmz*Hillszmz*Hillssmz*Hillsemz-kmrgd(2)*(0.04*6*Mu1+0.2*15*Mu2+20*Mu3+15*Mu4+6*Mu5+Mu6)-kmz*kmrgd(2);
gz*kmrgd(2)*(Mu0+0.6*6*Mu1+0.3*15*Mu2+0.1*20*Mu3+0.05*15*Mu4+0.05*6*Mu5+0.05*Mu6)-kz*kmrgd(3)
gs*HillsIs*Hillsss*Hillssls*Hillsks*Hillsws-ks*kmrgd(4);
gmsl*Hillssmsl*Hillskmsl*Hillswmsl-kmrgd(5)*(0.04*6*Mu1+0.2*15*Mu2+20*Mu3+15*Mu4+6*Mu5+Mu6)-kmsl*kmrgd(5);
gsl*kmrgd(5)*(Mu0+0.6*6*Mu1+0.3*15*Mu2+0.1*20*Mu3+0.05*15*Mu4+0.05*6*Mu5+0.05*Mu6)-ksl*kmrgd(6);
gk*Hillsslk*Hillssk*Hillskk-kk*kmrgd(7);
ge*Hillsse*Hillssle*Hillsee-ke*kmrgd(8);
gw*Hillsww-kw*kmrgd(9)
] ;