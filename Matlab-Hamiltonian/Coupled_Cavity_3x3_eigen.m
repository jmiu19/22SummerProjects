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
Lx=5; %exciton loss

gain=linspace(0,11,111);
loss=linspace(0,-11,111);

%% Polariton gain and loss tuning
amax=401; 
EpGL=zeros(3,amax);
  C=5;%Coupling
  GG=0;%Loss & Gain
  LL=0 %Loss
  Xl=0.01; %exciton loss
  XC=1650;
  detuning = 7;
  Ecm=XC+detuning;
  delta=50; 
  Ecc=Ecm-delta;
  RabiC=10;
  Rabic=10;
  
  

for a= 1:1:amax;
    G=GG+0.1.*(a-1);
    %G=0;
    %L=0;
    L=LL;
   % Photonic mode
    %Ecm=Ectm*(1-(sind(0))^(2)/ntm)^(-1/2);
    syms x
        M=[
            Ecm+1i*G,	           C,      RabiC/2;	         
                     C,   Ecm-1i*G,      Rabic/2;    
               RabiC/2,	     Rabic/2,	      XC
           ];
        [U V] = eig(M);
        EpGL(:,a)=[V(1,1); V(2,2); V(3,3)];
        eigv1(:,a)=U(:,1);
        eigv2(:,a)=U(:,2);
        eigv3(:,a)=U(:,3);
                disp(a)  
        
end
%% Extraction eigenvectors
    RealEV1=real(eigv1);
    ImagEV1=imag(eigv1);
    RealEV11=RealEV1(1:1,:);
    RealEV12=RealEV1(2:2,:);
    RealEV13=RealEV1(3:3,:);
    ImagEV11=ImagEV1(1:1,:);
    ImagEV12=ImagEV1(2:2,:);
    ImagEV13=ImagEV1(3:3,:);
    RealEV2=real(eigv2);
    ImagEV2=imag(eigv2);
    RealEV21=RealEV2(1:1,:);
    RealEV22=RealEV2(2:2,:);
    RealEV23=RealEV2(3:3,:);
    ImagEV21=ImagEV2(1:1,:);
    ImagEV22=ImagEV2(2:2,:);
    ImagEV23=ImagEV2(3:3,:);
    RealEV3=real(eigv3);
    ImagEV3=imag(eigv3);
    RealEV31=RealEV3(1:1,:);
    RealEV32=RealEV3(2:2,:);
    RealEV33=RealEV3(3:3,:);
    ImagEV31=ImagEV3(1:1,:);
    ImagEV32=ImagEV3(2:2,:);
    ImagEV33=ImagEV3(3:3,:);
%% Extraction eigenvalues
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
GL = GG-LL;
Gain= GL:0.1:cont; % contrast
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

%% Eigenvector Plot for vector 1 Real
figure(8)
plot(Gain,RealEV11,'r.','LineWidth',0.3);
hold on;
plot(Gain,RealEV12,'g.','LineWidth',0.3);
hold on;
plot(Gain,RealEV13,'b.','LineWidth',0.3);
hold on;
%% Eigenvector Plot for vector 1 Imag
figure(9)
plot(Gain,ImagEV11,'r.','LineWidth',0.3);
hold on;
plot(Gain,ImagEV12,'g.','LineWidth',0.3);
hold on;
plot(Gain,ImagEV13,'b.','LineWidth',0.3);
hold on;
