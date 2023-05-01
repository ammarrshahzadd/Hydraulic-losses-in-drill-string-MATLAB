clc
clear all
%Defining constants
q = [100 200 300 400 500 600 700 800 900 1000];
IDdp=3.5;
IDdc=4;
IDhole=6.5;
IDcasing=7.5;
ODdc=5.5;
ODdp=4.5;
den= 10.5;
Theta300=36;
Theta600=60;
E = 0.00065;
L_DP=6450;
L_DC=550;
L_DCA=550;
L_DPA1=950;
L_DPA2=5500;
Cd=0.95;
DPsurface=50;
ROP=40;
por=0.2;
Dcut=0.025;
gammarock=2.7;
gammamud = 10.5;
Vslip = [1 1 1 1 1 1 1 1 1 1];

%Plotting critical reynolds number plot
HN = [1e3 2e3 3e3 4e3 5e3 6e3 7e3 8e3 9e3 1e4 2e4 5e4 1e5 2e5 3e5 1e6 5e6 1e7];
CRe = [2.2e3 2.5e3 2.6e3 2.75e3 2.9e3 3e3 3.1e3 3.2e3 3.3e3 3.4e3 4e3 5.5e3 6.9e3 8.9e3 1e4 1.6e4 2.85e4 3.7e4];
loglog(HN,CRe)
ylim([1e3 1e5])
grid on

%Bingham plastic Model
for counter=1:1:10
    n=1;
%Mean velocity
%Drill Pipe
v_DP(counter)=q(counter)/(2.448*IDdp^2);

%Drill collar
v_DC(counter)=q(counter)/(2.448*IDdc^2);

%annulus
%section 1
%Drill collar annulus
v_DCA(counter)=q(counter)/(2.448*(IDhole^(2)-ODdc^(2)));
%section 2
%Drill pipe annulus
v_DPA1(counter)=q(counter)/(2.448*(IDhole^(2)-ODdp^(2)));
%section 3
%casing annulus
v_DPA2(counter)=q(counter)/(2.448*(IDcasing^(2)-ODdp^(2)));

%flow behavior parameter
PlasticVis=Theta600-Theta300;
Tao=Theta300-PlasticVis;

%Turbulence criteria
%Drill pipe
Nhep=(37100*den*Tao*IDdp^(2))/PlasticVis^(2);
Nrecp = interp1(HN,CRe,Nhep);
Nrep(counter)=(928*den*v_DP(counter)*IDdp)/PlasticVis;
fp(counter) = ((0.25)/(log10(((E/IDdp)/3.715)+(6.943/Nrep(counter))^0.9))^2)/4;

if Nrep(counter)>Nrecp
    dPdivdL_p(counter)=(fp(counter)*den*v_DP(counter)^(2))/(25.8*IDdp);
else 
    dPdivdL_p(counter)=((PlasticVis*v_DP(counter))/(1500*IDdp^2))+(Tao/(225*IDdp));
end

dp_DP(counter)=dPdivdL_p(counter)*L_DP;

%Drill collar
Nhec=(37100*den*Tao*IDdc^(2))/PlasticVis^(2);
Nrecc = interp1(HN,CRe,Nhec);
Nrec(counter)=(928*den*v_DC(counter)*IDdc)/PlasticVis;
fc(counter) = ((0.25)/(log10(((E/IDdc)/3.715)+(6.943/Nrec(counter))^0.9))^2)/4;

if Nrec(counter)>Nrecc
    dPdivdL_c(counter)=(fc(counter)*den*v_DC(counter)^(2))/(25.8*IDdc);
else 
    dPdivdL_c(counter)=((PlasticVis*v_DC(counter))/(1500*IDdc^2))+(Tao/(225*IDdc));
end

dp_DC(counter)=dPdivdL_c(counter)*L_DC;

%Annulus
%section 1
%Drill collar annulus
Nheca=(24700*den*Tao*(IDhole-ODdc)^2)/PlasticVis^(2);
Nrecca = interp1(HN,CRe,Nheca);
Nreca(counter)=(757*den*v_DCA(counter)*(IDhole-ODdc))/PlasticVis;
fca(counter) = ((0.25)/(log10(((E/(IDhole-ODdc))/3.715)+(6.943/Nreca(counter))^0.9))^2)/4;

if Nreca(counter)>Nrecca
    dPdivdL_ca(counter)=(fca(counter)*den*v_DCA(counter)^(2))/(21.1*(IDhole-ODdc));
else 
    dPdivdL_ca(counter)=((PlasticVis*v_DCA(counter))/(1000*(IDhole-ODdc)^2))+(Tao/(200*(IDhole-ODdc)));
end

dp_DCA(counter)=dPdivdL_ca(counter)*L_DCA;

%section 2
%Drill pipe annulus
Nhepa=(24700*den*Tao*(IDhole-ODdp)^2)/PlasticVis^(2);
Nrecpa = interp1(HN,CRe,Nhepa);
Nrepa(counter)=(757*den*v_DPA1(counter)*(IDhole-ODdp))/PlasticVis;

fpa(counter) = ((0.25)/(log10(((E/(IDhole-ODdp))/3.715)+(6.943/Nrepa(counter))^0.9))^2)/4;

if Nrepa(counter)>Nrecpa
    dPdivdL_pa(counter)=(fpa(counter)*den*v_DPA1(counter)^(2))/(21.1*(IDhole-ODdp));
else 
    dPdivdL_pa(counter)=((PlasticVis*v_DPA1(counter))/(1000*(IDhole-ODdp)^2))+(Tao/(200*(IDhole-ODdp)));
end

dp_DPA1(counter)=dPdivdL_pa(counter)*L_DPA1;

%section 3
%Casing annulus
Nhea=(24700*den*Tao*(IDcasing-ODdp)^2)/PlasticVis^(2);
Nreca = interp1(HN,CRe,Nhea);
Nrea(counter)=(757*den*v_DPA2(counter)*(IDcasing-ODdp))/PlasticVis;

fa(counter) = ((0.25)/(log10(((E/(IDcasing-ODdp))/3.715)+(6.943/Nrea(counter))^0.9))^2)/4;

if Nrea(counter)>Nreca
    dPdivdL_a(counter)=(fa(counter)*den*v_DPA2(counter)^(2))/(21.1*(IDcasing-ODdp));
else 
    dPdivdL_a(counter)=((PlasticVis*v_DPA2(counter))/(1000*(IDcasing-ODdp)^2))+(Tao/(200*(IDcasing-ODdp)));
end

dp_DPA2(counter)=dPdivdL_a(counter)*L_DPA2;

%Pressure drop at Drill bit
%Three cone roller bit
D1=15/32;
D2=12/32;
D3=15/32;

At=(pi/4)*(D1^(2)+D2^(2)+D3^(2));

dPbit(counter)=((8.311e-5)*den*q(counter)^2)/(Cd^(2)*At^2);

%Total Pressure drop
DPtotal(counter)=dp_DP(counter)+dp_DC(counter)+dp_DCA(counter)+dp_DPA1(counter)+dp_DPA2(counter)+dPbit(counter)+DPsurface;


%Cuttings Transport
A_bit=(3.1428/4)*(IDhole/12)^2;
VolCut=(pi()/6)*(Dcut^3);

w = VolCut*gammarock*8.33;
fbuoy = VolCut*gammamud;
ApVis = PlasticVis + (5*Tao*Dcut)/v_DCA(counter);

fdrag=zeros(1,10);
while(w ~= fbuoy + fdrag(counter) && n<1000)
NreParticle(counter) = 928*gammamud*Vslip(counter)*Dcut/ApVis;
DragCoeff(counter) = (24/NreParticle(counter))+(3/NreParticle(counter)^0.5)+0.34;
Vslip(counter) = 1.89*((Dcut/DragCoeff(counter))*((gammarock*8.33-gammamud)/(gammamud)))^0.5;
fdrag(counter) = DragCoeff(counter)*(0.5*gammamud*Vslip(counter)^2)*(pi/4)*Dcut^2;
n=n+1;
end
v_Cut(counter) = v_DCA(counter)-Vslip(counter);

A_bit=At;
Q_rock=ROP*A_bit*(1-por);
A_annulus=((pi/4)*(IDhole^2-ODdc^2))+(pi/4)*((IDhole^2-ODdp^2))+(pi/4)*((IDcasing^2-ODdp^2));
f_Cutt(counter)=Q_rock/(A_annulus*v_Cut(counter));
rho_mud(counter)=(f_Cutt(counter)*gammarock)+(1-f_Cutt(counter))*gammamud;

end
%Plotting
f1=figure;
f2=figure;
f3=figure;
f4=figure;
f5=figure;

%q vs v_DCut
figure(f1);
plot(q,v_Cut)
xlim([0 1000])
legend({'velcutt'},'Location','northwest')
grid on
%q Vs Pressure loss
figure(f2);
plot(q,dp_DPA2)
xlim([0 1000])
hold on
plot(q,dp_DP)
plot(q,dp_DC)
plot(q,dPbit)
plot(q,dp_DCA)
plot(q,dp_DPA1)
legend({'Pressure Loss Casing Annulus','Pressure Loss Pipe','Pressure Loss Drill Collar','Pressure Loss Bit','Pressure Loss Drill Collar Annulus','Pressure Loss Drill Pipe Annulus'},'Location','northwest')
hold off
grid on
%q vs Total frictional pressure loss
figure(f3);
plot(q,DPtotal)
xlim([0 1000])
legend({'DPtotal'},'Location','northwest')
grid on
%q vs f_cutt
figure(f4);
plot(q,f_Cutt)
xlim([0 1000])
legend({'f_Cutt'},'Location','northwest')
grid on
%q vs rho_mud
figure(f5);
plot(q,rho_mud)
xlim([0 1000])
legend({'rho_mud'},'Location','northwest')
grid on