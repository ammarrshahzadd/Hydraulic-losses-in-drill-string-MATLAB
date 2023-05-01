 
%Newtonian Model
%Mean velocity
%Drill Pipe
q=400;
IDdp=3.5;
IDdc=4;
IDhole=6.5;
IDcasing=7.5;
ODdc=5.5;
ODdp=4.5;
den=10.5;
Theta300=36;
f=0.025;


v=q/(2.448*IDdp^(2));
%Drill collar
v=q/(2.448*IDdc^(2));
%annulus
%section 1
%Drill collar annulus
Vis=q/(2.448*(IDhole^(2)-ODdc^(2)));
%section 2
%Drill pipe annulus
v=q/(2.448*(IDhole^(2)-ODdp^(2)));
%section 3
%casing annulus
v=q/(2.448*(IDcasing^(2)-ODdp^(2)));
%flow behavior parameter
Vis=Theta300;
%Turbulence criteria
%Drill pipe
Nrec=2100;    
Nre=(928*den*v*IDdp)/Vis;
if Nre>Nrec
    dPdivdL=(f*den*v^(2))/(25.8*IDdp);
else 
   dPdivdL=(vis*v)/(1500*IDdp^2);
end    
%Drill collar
Nrec=2100;    
Nre=(928*den*v*IDdc)/Vis;
if Nre>Nrec
    dPdivdL=(f*den*v^(2))/(25.8*IDdc);
else 
    dPdivdL=(vis*v)/(1500*IDdc^2);
end    
%Annulus
%section 1
%Drill collar annulus
Nrec=2100;
if Nre>Nrec
    dPdivdL=(f*den*v)/(21.1*(IDhole-ODdc));
else 
    dPdivdL=(vis*v)/(1000*(IDhole-ODdc)^2);
end  
Nre=(757*den*v*(IDhole-ODdc))/Vis;
%section 2
%Drill pipe annulus
Nrec=2100;
if Nre>Nrec
    dPdivdL=(f*den*v)/(21.1*(IDhole-ODdp));
else 
  dPdivdL=(vis*v)/(1000*(IDhole-ODdc)^2);
end  
Nre=(757*den*v*(IDhole-ODdp))/vis;
%section 3
%Casing annulus
Nrec=2100;
if Nre>Nrec
   dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdp));
else 
    dPdivdL=(vis*v)/(1000*(IDcasing-ODdp)^2);
end  
Nre=(757*den*v*(IDcasing-ODdp))/vis;
%Laminar Flow frictional pressure loss
%Drill Pipe
dPdivdL=(vis*v)/(1500*IDdp^2);
%Drill collar
dPdivdL=(vis*v)/(1500*IDdc^2);
%Annulus
%section 1
%Drill collar annulus
dPdivdL=(vis*v)/(1000*(IDcasing-ODdc)^2);
%section 2
%Drill pipe annulus
dPdivdL=(vis*v)/(1000*(IDcasing-ODdp)^2);
%section 3
%Casing annulus
dPdivdL=(vis*v)/(1000*(IDcasing-ODdp)^2);
%Turbulent Flow frictional Pressure loss
%Drill Pipe
dPdivdL=(f*den*v^(2))/(25.8*IDdp);
%Drill collar
dPdivdL=(f*den*v^(2))/(25.8*IDdc);
%Annulus
%section 1
%Drill collar annulus
dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdc));
%section 2
%Drill Pipe annulus
dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdp));
%section 3
%Casing annulus
dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdp));

%Bingham plastic Model
%Mean velocity
%Drill Pipe
V=q/(2.448*IDdp^2);
%Drill collar
V=q/(2.448*IDdc^2);
%annulus
%section 1
%Drill collar annulus
v=q/(2.448*(IDhole^(2)-ODdc^(2)));
%section 2
%Drill pipe annulus
v=q/(2.448*(IDhole^(2)-ODdp^(2)));
%section 3
%casing annulus
v=q/(2.448*(IDcasing^(2)-ODdp^(2)));
%flow behavior parameter
PlasticVis=Theta600-Theta300;
Tao=Theta300-PlasticVis;
%Turbulence criteria
%Drill pipe
Nhe=(37100*den*tao*IDdp^(2))/PlasticVis^(2);

if Nrec==2100
    apparentVis=PlasticVis+((6.66*tao*IDdp)/v);
    PlasticVis=apparentVis;
elseif Nrec>Nre
    PlasticVis=Theta600-Theta300;
end
    
Nre=(928*den*v*IDdp)/PlasticVis;
%Drill collar
Nhe=(37100*den*tao*IDdc^(2))/PlasticVis^(2);

if Nrec==2100
    apparentVis=PlasticVis+((6.66*tao*IDdc)/v);
    PlasticVis=apparentVis;
elseif Nrec>Nre
    PlasticVis=Theta600-Theta300;
end
    
Nre=(928*den*v*IDdc)/PlasticVis;
%Annulus
%section 1
%Drill collar annulus
Nhe=(24700*den*tao*(IDcasing-ODdc)^2)/PlasticVis^(2);

if Nrec==2100
    apparentVis=PlasticVis+((5*tao*(IDcasing-ODdc))/v);
    
    %Not sure how to use if else statment here for now
    
Nre=(757*den*tao*(IDcasing-ODdc))/PlasticVis;
%section 2
%Drill pipe annulus
Nhe=(24700*den*tao*(IDcasing-ODdp)^2)/PlasticVis^(2);

if Nrec==2100
    apparentVis=PlasticVis+((5*tao*(IDcasing-ODdp))/v);
    
     %Not sure how to use if else statment here for now
    
Nre=(757*den*tao*(IDcasing-ODdp))/PlasticVis;
%section 3
%Casing annulus
Nhe=(24700*den*tao*(IDcasing-ODdp)^2)/PlasticVis^(2);

if Nrec==2100
    apparentVis=PlasticVis+((5*tao*(IDcasing-ODdp))/v);
    
     %Not sure how to use if else statment here for now
    
Nre=(757*den*tao*(IDcasing-ODcasing))/PlasticVis;
%Laminar Flow frictional pressure loss
%Drill Pipe
dPdivdL=((PlasticVis*v)/(1500*IDdp^2))+(tao/(225*IDdp));
%Drill collar
dPdivdL=((PlasticVis*v)/(1500*IDdc^2))+(tao/(225*IDdc));
%Annulus
%section 1
%Drill collar annulus
dPdivdL=((PlasticVis*v)/(1000*(IDcasing-ODdc)^2))+(tao/(200*(IDcasing-ODdc)));
%section 2
%Drill pipe annulus
dPdivdL=((PlasticVis*v)/(1000*(IDcasing-ODdp)^2))+(tao/(200*(IDcasing-ODdp)));
%section 3
%Casing annulus
dPdivdL=((PlasticVis*v)/(1000*(IDcasing-ODdp)^2))+(tao/(200*(IDcasing-ODdp)));
%Turbulent Flow frictional Pressure loss
%Drill Pipe
dPdivdL=(f*den*v^(2))/(25.8*IDdp);
%Drill collar
dPdivdL=(f*den*v^(2))/(25.8*IDdc);
%Annulus
%section 1
%Drill collar annulus
dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdc));
%section 2
%Drill Pipe annulus
dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdp));
%section 3
%Casing annulus
dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdp));

%Power law model
%Mean velocity
%Drill Pipe
v=q/(2.448*IDdp^2);
%Drill collar
v=q/(2.448*IDdc^2);
%annulus
%section 1
%Drill collar annulus
v=q/(2.448*(IDcasing^(2)-ODdc^(2)));
%section 2
%Drill pipe annulus
v=q/(2.448*(IDcasing^(2)-ODdp^(2)));
%section 3
%casing annulus
v=q/(2.448*(IDcasing^(2)-ODdp^(2)));
%flow behavior parameter
n=3.32*LOG(Theta600/Theta300);
K=(510*Theta300)/(511^n);
%Turbulence criteria
%Drill pipe

Nre=((89100*den*v^(2-n))/K)*((0.0416*IDdp)/(3+(1/n)))^n;

if Nre>Nrec
    dPdivdL=(f*den*v^(2))/(25.8*IDdp);
else 
    dPdivdL=(K*v^(n)*((3+(1/n))/0.0416)^n)/(144000*IDdp^(1+n));
end
%Drill collar

Nre=((89100*den*v^(2-n))/K)*((0.0416*IDdc)/(3+(1/n)))^n;

if Nre>Nrec
    dPdivdL=(f*den*v^(2))/(25.8*IDdc);
else 
    dPdivdL=(K*v^(n)*((3+(1/n))/0.0416)^n)/(144000*IDdc^(1+n));
end
%Annulus
%section 1
%Drill collar annulus

Nre=((109000*den*v^(2-n))/K)*((0.0208*(IDcasing-IDdc))/(2+(1/n)))^n;
if Nre>Nrec
    dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdc));
else 
    dPdivdL=(K*v^(n)*((2+(1/n))/0.0208)^n)/(144000*(IDcasing-IDdc)^(1+n));
end
%section 2
%Drill pipe annulus

Nre=((109000*den*v^(2-n))/K)*((0.0208*(IDcasing-IDdp))/(2+(1/n)))^n;
if Nre>Nrec
    dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdp));
else 
    dPdivdL=(K*v^(n)*((2+(1/n))/0.0208)^n)/(144000*(IDcasing-IDdp)^(1+n));
end
%section 3
%Casing annulus

Nre=((109000*den*v^(2-n))/K)*((0.0208*(IDcasing-ODdp))/(2+(1/n)))^n;
if Nre>Nrec
    dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdp));
else 
    dPdivdL=(K*v^(n)*((2+(1/n))/0.0208)^n)/(144000*(IDcasing-ODdp)^(1+n));
end
%Laminar Flow frictional pressure loss
%Drill Pipe
dPdivdL=(K*v^(n)*((3+(1/n))/0.0416)^n)/(144000*IDdp^(1+n));
%Drill collar
dPdivdL=(K*v^(n)*((3+(1/n))/0.0416)^n)/(144000*IDdc^(1+n));
%Annulus
%section 1
%Drill collar annulus
dPdivdL=(K*v^(n)*((2+(1/n))/0.0208)^n)/(144000*(IDcasing-IDdc)^(1+n));
%section 2
%Drill pipe annulus
dPdivdL=(K*v^(n)*((2+(1/n))/0.0208)^n)/(144000*(IDcasing-IDdp)^(1+n));
%section 3
%Casing annulus
dPdivdL=(K*v^(n)*((2+(1/n))/0.0208)^n)/(144000*(IDcasing-ODdp)^(1+n));
%Turbulent Flow frictional Pressure loss
%Drill Pipe
dPdivdL=(f*den*v^(2))/(25.8*IDdp);
%Drill collar
dPdivdL=(f*den*v^(2))/(25.8*IDdc);
%Annulus
%section 1
%Drill collar annulus
dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdc));
%section 2
%Drill Pipe annulus
dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdp));
%section 3
%Casing annulus
dPdivdL=(f*den*v)/(21.1*(IDcasing-ODdp));

%Still have to write the codes for drill bit 
%Need the graph values for Nrec plots for Bingham plastic and Power law 








