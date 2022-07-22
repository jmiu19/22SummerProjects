clear all;
close all;

  %% Costants
h=4.1357e-15;
c=299792458;
thmin= -68;
thmax= 68;
RabiA=40;
RabiB=40;
RabiC=40;
nte=7.3;
ntm=8.9;
R=8.5e-07;
%XA=3429;
%XB=3414;
XC=1653;
C=10;%Coupling experiment
gain=linspace(0,11,111);
loss=linspace(0,-11,111);
Ecm=1653;
Ecc=1630;
RabiC=20;

%% Polariton gain and loss tuning
amax=401; 
EpGL=zeros(3,amax);
  C=4;%Coupling
  G=10;%Loss & Gain
  L=10 %Loss
  Xl=0.01; %exciton loss
  Ecm=1690;
  delta=50; 
  Ecc=Ecm-delta;
  RabiC=40;
  XC=1620

for a= 1:1:amax;
    G=5+0.1.*(a-1);
    L=5;
   % Photonic mode
    %Ecm=Ectm*(1-(sind(0))^(2)/ntm)^(-1/2);
    syms x
        M=[
            Ecm+1i*G-x,	           C,      RabiC/2;	         
                     C,   Ecm-1i*L-x,      RabiC/2;    
               RabiC/2,	     RabiC/2,	      XC-x
           ];
        EpGL(:,a)=solve(det(M),x); 
end
%% Extraction
    RealEP=real(EpGL);
    ImagEP=imag(EpGL);
    RealGP=RealEP(1:1,:);
    RealLP=RealEP(2:2,:);
    RealEx=RealEP(3:3,:);
    %RealUP2=RealEP(4:4,:);
    ImagGP=ImagEP(1:1,:);
    ImagLP=ImagEP(2:2,:);
    ImagEx=ImagEP(3:3,:);
    %ImagUP2=ImagEP(4:4,:);

%% Real part Plot
figure(6)
cont = G-L;
Gain= 0:0.1:cont; % contrast
plot(Gain,RealGP,'r.','LineWidth',0.3);
hold on;
plot(Gain,RealLP,'g.','LineWidth',0.3);
hold on;
plot(Gain,RealEx,'b.','LineWidth',0.3);
hold on;
%plot(Gain,RealUP2,'b.','LineWidth',0.3);
%hold on;
%% Imaginaryt part Plot
figure(7)
plot(Gain,ImagGP,'r.','LineWidth',0.3);
hold on;
plot(Gain,ImagLP,'g.','LineWidth',0.3);
hold on;
plot(Gain,ImagEx,'b.','LineWidth',0.3);
hold on;
%plot(Gain,ImagUP2,'r.','LineWidth',0.3);