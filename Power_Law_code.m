clc
clear all

%Defining constants
q = [7.67 10 50 100 200 300 400 500 600 700 800];
IDdp=3.5;
IDdc=3.5;
IDhole=5.5;
IDcasing=6.5;
ODdc=4.5;
ODdp=4;
den=11;
Theta300=41;
Theta600=58;
E = 0.00065;
Ldp=2400;
Ldc=600;
Ldca=600;
Ldpa=600;
Lda=1800;
Cd=0.95;
DPsurface=30;
ROP=40;
por=0.18;
Dcut=0.025;
gammarock=2.7;
gammamud = 10.5;
PlasticVis=24;
Tao=12;
Vslip = [1 1 1 1 1 1 1 1 1 1 1];

%Power law model
for counter=1:1:11
    m=1;
%Mean velocity
%Drill Pipe
vp(counter)=q(counter)/(2.448*IDdp^2);
%Drill collar
vc(counter)=q(counter)/(2.448*IDdc^2);
%annulus
%section 1
%Drill collar annulus
vca(counter)=q(counter)/(2.448*(IDhole^(2)-ODdc^(2)));
%section 2
%Drill pipe annulus
vpa(counter)=q(counter)/(2.448*(IDhole^(2)-ODdp^(2)));
%section 3
%casing annulus
va(counter)=q(counter)/(2.448*(IDcasing^(2)-ODdp^(2)));

%flow behavior parameter
n=3.32*log10(Theta600/Theta300);
K=(510*Theta300)/(511^n);

%Turbulence criteria
%Drill pipe

Nrep(counter)=((89100*den*vp(counter)^(2-n))/K)*((0.0416*IDdp)/(3+(1/n)))^n;


if n>0.4
    Nrecp=2100;
else 
    Nrecp=407.04*n^(-0.906);
end

if Nrep(counter)>Nrecp
fp(counter)=0.00001;
a(counter) = (((4)/(n^0.75))*log10(Nrep(counter)*(fp(counter)^(1-n/2)))-0.395/(n^1.2))-(1/fp(counter))^0.5;

while (round(a(counter),1)~=0)
    fp(counter) = fp(counter)+0.00001;
a(counter) = (((4)/(n^0.75))*log10(Nrep(counter)*(fp(counter)^(1-n/2)))-0.395/(n^1.2))-(1/fp(counter))^0.5;
end
    dPdivdL_p(counter)=(fp(counter)*den*vp(counter)^(2))/(25.8*IDdp);
else 
    dPdivdL_p(counter)=(K*vp(counter)^(n)*((3+(1/n))/0.0416)^n)/(144000*IDdp^(1+n));
end
dp_DP(counter)=dPdivdL_p(counter)*Ldp;


%Drill collar
Nrec(counter)=((89100*den*vc(counter)^(2-n))/K)*((0.0416*IDdc)/(3+(1/n)))^n;

if n>0.4
    Nrecc=2100;
else 
    Nrecc=407.04*n^(-0.906);
end

if Nrec(counter)>Nrecc

fc(counter) = 0.00001;
b(counter) = fc(counter)*(((4)/(n^0.75))*log10(Nrec(counter)*fc(counter)^(1-n/2))-0.395/(n^1.2))^2-sqrt(1/fc(counter));
while(round(b(counter),2)~=0)
b(counter) = fc(counter)*(((4)/(n^0.75))*log10(Nrec(counter)*fc(counter)^(1-n/2))-0.395/(n^1.2))^2-sqrt(1/fc(counter));
fc(counter)=fc(counter)+0.00001;
end
    dPdivdL_c(counter)=(fc(counter)*den*vc(counter)^(2))/(25.8*IDdc);
    
else 
    dPdivdL_c(counter)=(K*vc(counter)^(n)*((3+(1/n))/0.0416)^n)/(144000*IDdc^(1+n));
end
dp_DC(counter)=dPdivdL_c(counter)*Ldc;
%Annulus
%section 1
%Drill collar annulus

Nreca(counter)=((109000*den*vca(counter)^(2-n))/K)*((0.0208*(IDhole-ODdc))/(2+(1/n)))^n;

if n>0.4
    Nrecca=2100;
else 
    Nrecca=407.04*n^(-0.906);
end

if Nreca(counter)>Nrecca
fca(counter) = 0.00001;
c(counter) = fca(counter)*(((4)/(n^0.75))*log10(Nreca(counter)*fca(counter)^(1-n/2))-0.395/(n^1.2))^2-sqrt(1/fca(counter));
while(round(c(counter),2)~=0)
c(counter) = fca(counter)*(((4)/(n^0.75))*log10(Nreca(counter)*fca(counter)^(1-n/2))-0.395/(n^1.2))^2-sqrt(1/fca(counter));
fca(counter)=fca(counter)+0.00001;
end
    dPdivdL_ca(counter)=(fca(counter)*den*vca(counter))/(21.1*(IDhole-ODdc));
else 
    dPdivdL_ca(counter)=(K*vca(counter)^(n)*((2+(1/n))/0.0208)^n)/(144000*(IDhole-ODdc)^(1+n));
end
dp_DCA(counter)=dPdivdL_ca(counter)*Ldca;
%section 2
%Drill pipe annulus

Nrepa(counter)=((109000*den*vpa(counter)^(2-n))/K)*((0.0208*(IDhole-ODdp))/(2+(1/n)))^n;

if n>0.4
    Nrecpa=2100;
else 
    Nrecpa=407.04*n^(-0.906);
end

if Nrepa(counter)>Nrecpa
fpa(counter) = 0.00001;
d(counter) = fpa(counter)*(((4)/(n^0.75))*log10(Nrepa(counter)*fpa(counter)^(1-n/2))-0.395/(n^1.2))^2-sqrt(1/fpa(counter));
while(round(d,2)~=0)
d(counter) = fpa(counter)*(((4)/(n^0.75))*log10(Nrepa(counter)*fpa(counter)^(1-n/2))-0.395/(n^1.2))^2-sqrt(1/fpa(counter));
fpa(counter)=fpa(counter)+0.00001;
end
    dPdivdL_pa(counter)=(fpa(counter)*den*vpa(counter))/(21.1*(IDhole-ODdp));
else 
    dPdivdL_pa(counter)=(K*vpa(counter)^(n)*((2+(1/n))/0.0208)^n)/(144000*(IDhole-ODdp)^(1+n));
end
dp_DPA1(counter)=dPdivdL_pa(counter)*Ldpa;

%section 3
%Casing annulus

Nrea(counter)=((109000*den*va(counter)^(2-n))/K)*((0.0208*(IDcasing-ODdp))/(2+(1/n)))^n;

if n>0.4
    Nreca=2100;
else 
    Nreca=407.04*n^(-0.906);
end

if Nrea(counter)>Nreca
fa(counter) = 0.00001;
e(counter) = fa(counter)*(((4)/(n^0.75))*log10(Nrea(counter)*fa(counter)^(1-n/2))-0.395/(n^1.2))^2-sqrt(1/fa(counter));
while(round(e,2)~=0)
e(counter) = fa(counter)*(((4)/(n^0.75))*log10(Nrea(counter)*fa(counter)^(1-n/2))-0.395/(n^1.2))^2-sqrt(1/fa(counter));
fa(counter)=fa(counter)+0.00001;
end
    dPdivdL_a(counter)=(fa(counter)*den*va(counter))/(21.1*(IDcasing-ODdp));
else 
    dPdivdL_a(counter)=(K*va(counter)^(n)*((2+(1/n))/0.0208)^n)/(144000*(IDcasing-ODdp)^(1+n));
end
dp_DPA2(counter)=dPdivdL_a(counter)*Lda;

%Pressure drop at Drill bit
%Three cone roller bit
D1=12/32;
D2=12/32;
D3=12/32;

At=(pi/4)*(D1^(2)+D2^(2)+D3^(2));

dPbit(counter)=((8.311e-5)*den*q(counter)^2)/(Cd^(2)*At^2);

%Total Pressure drop
DPtotal(counter)=dp_DP(counter)+dp_DC(counter)+dp_DCA(counter)+dp_DPA1(counter)+dp_DPA2(counter)+dPbit(counter)+DPsurface;

%Cuttings Transport
A_bit=(3.1428/4)*(IDhole/12)^2;
VolCut=(pi()/6)*(Dcut^3);

w = VolCut*gammarock*8.33;
fbuoy = VolCut*gammamud;
ApVis = PlasticVis + (5*Tao*Dcut)/vca(counter);

fdrag=zeros(1,12);
while(w ~= fbuoy + fdrag(counter) && m<1000)
NreParticle(counter) = 928*gammamud*Vslip(counter)*Dcut/ApVis;
DragCoeff(counter) = (24/NreParticle(counter))+(3/NreParticle(counter)^0.5)+0.34;
Vslip(counter) = 1.89*((Dcut/DragCoeff(counter))*((gammarock*8.33-gammamud)/(gammamud)))^0.5;
fdrag(counter) = DragCoeff(counter)*(0.5*gammamud*Vslip(counter)^2)*(pi/4)*Dcut^2;
m=m+1;
end
vCut(counter) = vca(counter)-Vslip(counter);
A_bit=At;
Q_rock=ROP*A_bit*(1-por);
A_annulus=((pi/4)*(IDhole^2-ODdc^2))+(pi/4)*((IDhole^2-ODdp^2))+(pi/4)*((IDcasing^2-ODdp^2));
f_Cutt(counter)=Q_rock/(A_annulus*vCut(counter));
rho_mud(counter)=(f_Cutt(counter)*gammarock)+(1-f_Cutt(counter))*gammamud;

end


%Plotting
f1=figure;
f2=figure;
f3=figure;
f4=figure;
f5=figure;


%q vs Vcut
figure(f1);
plot(q,vCut)
xlim([0 800])
legend({'velcutt'},'Location','northwest')
grid on

%q Vs Pressure loss
figure(f2);
plot(q,dp_DPA2)
xlim([0 1000])
grid on
hold on
plot(q,dp_DP)
plot(q,dp_DC)
plot(q,dPbit)
plot(q,dp_DCA)
plot(q,dp_DPA1)
legend({'Pressure Loss Casing Annulus','Pressure Loss Pipe','Pressure Loss Drill Collar','Pressure Loss Bit','Pressure Loss Drill Collar Annulus','Pressure Loss Drill Pipe Annulus'},'Location','northwest')
hold off

%q vs Total frictional pressure loss
figure(f3);
plot(q,DPtotal)
xlim([0 800])
legend({'DPtotal'},'Location','northwest')
grid on
%q vs f_cutt
figure(f4);
plot(q,f_Cutt)
xlim([0 800])
legend({'f_Cutt'},'Location','northwest')
grid on
%q vs rho_mud
figure(f5);
plot(q,rho_mud)
xlim([0 800])
legend({'rho_mud'},'Location','northwest')
grid on

