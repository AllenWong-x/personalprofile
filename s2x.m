clear;
clc;
A= xlsread('Sdata.xlsx','','B1:b201');
B= xlsread('Sdata.xlsx','','C1:C201');
C= xlsread('Sdata.xlsx','','D1:D201');
D= xlsread('Sdata.xlsx','','E1:E201');
E= xlsread('Sdata.xlsx','','F1:F201');
F= xlsread('Sdata.xlsx','','G1:G201');
G= xlsread('Sdata.xlsx','','H1:H201');
H= xlsread('Sdata.xlsx','','I1:I201');
T= xlsread('Sdata.xlsx','','A1:A201');
for i = 1:201
    S11(i) = A(i).*cos(B(i)) + j*A(i).*sin(B(i));
    S21(i) = C(i).*cos(D(i)) + j*C(i).*sin(D(i));
    S12(i) = E(i).*cos(F(i)) + j*E(i).*sin(F(i));
    S22(i) = G(i).*cos(H(i)) + j*G(i).*sin(H(i));
    S(i,1)=T(i);
    S(i,2)=S11(i);
    S(i,3)=S21(i);
    S(i,4)=S12(i);
    S(i,5)=S22(i);
end

Z0=50;
for i=1:201
    Z11(i)=Z0.*[(1+S11(i)).*(1-S22(i))+S12(i).*S21(i)]./[(1-S11(i)).*(1-S22(i))-S12(i).*S21(i)];
    Z12(i)=Z0.*2.*S12(i)./[(1-S11(i)).*(1-S22(i))-S12(i).*S21(i)];
    Z21(i)=Z0.*2.*S21(i)./[(1-S11(i)).*(1-S22(i))-S12(i).*S21(i)];
    Z22(i)=Z0.*[(1-S11(i)).*(1+S22(i))+S12(i).*S21(i)]./[(1-S11(i)).*(1-S22(i))-S12(i).*S21(i)];
    
    Y11(i)=(1/Z0).*[(1-S11(i)).*(1+S22(i))+S12(i).*S21(i)]./[(1+S11(i)).*(1+S22(i))-S12(i).*S21(i)];
    Y12(i)=(1/Z0).*(-2).*S12(i)./[(1+S11(i)).*(1+S22(i))-S12(i)*S21(i)];
    Y21(i)=(1/Z0).*(-2).*S21(i)./[(1+S11(i)).*(1+S22(i))-S12(i).*S21(i)];
    Y22(i)=(1/Z0).*[(1+S11(i)).*(1-S22(i))+S12(i).*S21(i)]./[(1+S11(i)).*(1+S22(i))-S12(i).*S21(i)];
    
    Z(i,1)=T(i);
    Z(i,2)=Z11(i);
    Z(i,3)=Z12(i);
    Z(i,4)=Z21(i);
    Z(i,5)=Z22(i);
    
    Y(i,1)=T(i);
    Y(i,2)=Y11(i);
    Y(i,3)=Y12(i);
    Y(i,4)=Y21(i);
    Y(i,5)=Y22(i);
end   

for Cpgd=0:100
    for i=1:201
        w=0.15915494309189533576888376337251;
        Cgs(i)=2*Y12(i)./(-1j*w)-2*Cpgd;
        Cpg(i)=(Y11(i)+Y12(i))/(1j*w)-2*Cgs(i);
        Cds(i)=(Y22(i)+Y12(i))/(1j*w)-Cpg(i);
        Cpd(i)=Cpg(i);
        Cgd(i)=1/2*Cgs(i);
        Rpgd=1/(Y11(i)-1j*w*(Cpg(i)+Cgs(i)+Cgd(i)+Cpgd));
      W11(i)=1j*w*(Cpg(i)+2.5*Cgs(i)+Cpgd);
      W22(i)=1j*w*(Cpg(i)+Cds(i)+0.5*Cgs(i)+Cpgd);
      W21(i)=-1j*w*(0.5*Cgs(i)+Cpgd);
      W12(i)=-1j*w*(0.5*Cgs(i)+Cpgd);
    W(i,1)=W11(i);
    W(i,2)=W12(i);
    W(i,3)=W21(i);
    W(i,4)=W22(i);
    end
end 


for i=1:201 
    sigema(i)=W(i)-Y(i);
    if (min(sigema(i))<1*E-6)
    C(i,1)=Cpg(i);
    C(i,2)=Cpd(i);
    C(i,3)=Cds(i);
    C(i,4)=Cpgd;
    end
end
 
    
    figure(1)
    plot(T,real(Cpg)+imag(Cpg),'-r');
    figure(2)
    plot(T,real(Cgs)+imag(Cgs),'-b');
    figure(3)
    plot(T,real(Cds)+imag(Cds),'-g');
    figure(4)
    plot(T,real(sigema),'-black');
    
for i=1:201
     w=1.5915494309189533576888376337251e-10;
     Lg(i)=(imag(w*Z11(i))-imag(w*Z12(i)))/(w*w);
     Ld(i)=(imag(w*Z22(i))-imag(w*Z12(i)))/(w*w);
     Ls(i)=imag(w*Z12(i))/(w*w);
     
     Rg(i)=(real(w*Z11(i))-real(w*Z12(i)))/(w*w);
     Rd(i)=(real(w*Z22(i))-real(w*Z12(i)))/(w*w);
     Rs(i)=real(w*Z12(i))/(w*w);
     
    L(i,1)=T(i);
    L(i,2)=Lg(i);
    L(i,3)=Ld(i);
    L(i,4)=Ls(i);
    
    R(i,1)=T(i);
    R(i,2)=Rg(i);
    R(i,3)=Rd(i);
    R(i,4)=Rs(i);
    
end
