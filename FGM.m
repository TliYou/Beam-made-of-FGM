clc
clear all
syms x l z Ez Gz alpha_z
nel =input('Enter no of elements: ');
y = nel+1;
nnel=2;                                   % number of nodes per element
ndof=3;                                   % no. of degrees of freedom
nnode=nel+1;                              % total no. of nodes 
sdof=nnode*ndof;                          % total no. of dof  
n1=sdof;
L= 1/nel;
ele_l=1/nel;
Eb=390e9;
Gb=137e9;
alpha_b=6.9e-6;
Et=210e9;
Gt=80e9;
alpha_t=14e-6;
rhot=7850;
rhob=3950;
 
h=0.1;
n =input('Enter parameter of power law:  ');
width=0.1;
lim=0.1;
Ez=(Et-Eb)*(z/h+1/2)^n+Eb;              
Gz=(Gt-Gb)*(z/h+1/2)^n+Gb;
rhoz=(rhot-rhob)*(x/ele_l+1/2)^n+rhob;
alpha_z=(alpha_t-alpha_b)*(z/h+1/2)^n+alpha_b;
    I0=vpa(lim*int(rhoz,x,0,ele_l));
    I1=vpa(lim*int(x*rhoz,x,0,ele_l));
    I2=vpa(lim*int(x^2*rhoz,x,0,ele_l));
    A11=vpa(width*int(Ez,z,0,lim));
    B11=vpa(width*int((z*Ez),z,0,lim));
    D11=vpa(width*int((z^2*Ez),z,0,lim));
    A55=vpa(width*int(Gz,z,0,lim));
    %AT11=vpa(width*int((Ez*alpha_z),z,0,lim));
    %BT11=vpa(width*int((z*Ez*alpha_z),z,0,lim));
 
Mb=zeros(sdof,sdof);
Kb=zeros(sdof,sdof);
index=zeros(1,nnel*ndof);
%for n=0.5,cross section area = 0.001.
a = (B11*A55)/(A11*D11-B11^2);
b = A11*A55/(A11*D11-B11^2);
p = 1/(1+b*L^2/12);
 
% Mu
N1= [(1-x/l);6*a*p*((x/l)-1)*x;3*a*p*(x-l)*x;x/l;-6*a*p*((x/l)-1)*x;3*a*p*(x-l)*x];
 
N2= [(1-x/l) 6*a*p*((x/l)-1)*x 3*a*p*(x-l)*x x/l -6*a*p*((x/l)-1)*x 3*a*p*(x-l)*x];
 
M= I0*(N1*N2);
 
M1=vpa(int(M,x,0,1));
Mu=subs(M1,l,ele_l);
 
%Mw
N3=[0;
    p*(b*l^3+12*l-12*x+2*b*x^3-3*b*x^2*l)/l;
    x*p*(b*l^3+6*l-6*x+b*x^2*l-2*x*b*l^2)/l;
    0;
    -x*p*(-12+2*b*x^2-3*b*x*l)/l;
    x*p*(-6*l+b*x^2*l-x*b*l^2+6*x)/l];
 
N4=[0 p*(b*l^3+12*l-12*x+2*b*x^3-3*b*x^2*l)/l x*p*(b*l^3+6*l-6*x+b*x^2*l-2*x*b*l^2)/l 0 -x*p*(-12+2*b*x^2-3*b*x*l)/l x*p*(-6*l+b*x^2*l-x*b*l^2+6*x)/l];
 
M= I0*(N3*N4);
M2=vpa(int(M,x,0,1));
Mw=subs(M2,l,ele_l);
 
%Mphi
 
N5=[0;6*x*b*(-l+x)*p/l;(3*b*x^2*l+b*l^3+12*l-4*x*b*l^2-12*x)*p/l;0;-6*x*b*(-l+x)*p/l;x*(3*b*x*l-2*b*l^2+12)*p/l];
N6=[0 6*x*b*(-l+x)*p/l (3*b*x^2*l+b*l^3+12*l-4*x*b*l^2-12*x)*p/l 0 -6*x*b*(-l+x)*p/l x*(3*b*x*l-2*b*l^2+12)*p/l];
 
M= I2*(N5*N6);
M3=vpa(int(M,x,0,1));
Mphi=subs(M3,l,ele_l);
 
%Muphi
 
N7=(N1*N6)+(N5*N2);
M= -I1*N7;
M4=vpa(int(M,x,0,1));
Muphi=subs(M4,l,ele_l);
 
Mass_matrix = Mu + Mw + Mphi + Muphi;
 
for iel=1:nel
        index=sysdofel(iel,nnel,ndof);
        me=Mass_matrix;
% %         [k,z1,z2]=expKmtrx(lim);
%         k=DefIntKmtrx();
        Mb=assembly(Mb,me,index);
end
 
    
    
k1(1,1)=A11/L;
k1(1,3)=-B11/L;
k1(2,2)=A55*p/L;
k1(2,3)=A55*p/2;
k1(3,3)=D11/L+A55*p*L/4;
k1(3,6)=-D11/L+A55*p*L/4;
k1(1,4)=-k1(1,1);
k1(1,6)=-k1(1,3);
k1(2,5)=-k1(2,2);
k1(2,6)=k1(2,3);
k1(6,6)=k1(3,3);
k1(3,4)=-k1(1,3);
k1(3,5)=-k1(2,3);
k1(4,4)=k1(1,1);
k1(4,6)=k1(1,3);
k1(5,5)=k1(2,2);
k1(5,6)=-k1(2,3);   
for iel=1:nel
        index=sysdofel(iel,nnel,ndof);
        mk=k1;
% %         [k,z1,z2]=expKmtrx(lim);
%         k=DefIntKmtrx();
        Kb=assembly(Kb,k1,index);
end
 
lambda=eig(Kb,Mb);
natural_freq = sqrt(lambda);
%disp(natural_freq(1:10,1));
disp(natural_freq(3,1)/1e4);
disp(natural_freq(1,1)/1e3);
[V,D] = eig(Kb,Mb);
modeshapes=V;
v1=V(:,1);
v2=V(:,2);
v3=V(:,10);
 
first=zeros((n1/2)-1,1);
k=1;
for k=1:(n1/2)-1
    first(k) = first(k) + v1((2*k)-1);
    k=k+1;
end
 
sec=zeros((n1/2)-1,1);
k=1;
for k=1:(n1/2)-1
    sec(k) = sec(k) + v2((2*k)-1);
    k=k+1;
end
 
third=zeros((n1/2)-1,1);
k=1;
for k=1:(n1/2)-1
    third(k) = third(k) + v3((2*k)-1);
    k=k+1;
end
 
%MODE SHAPES for 1st three natural frequency displacement degree of freedom
subplot(3,2,1);
plot(first(:,1));
title('MODE SHAPE for 1st natural frequency(displacement)')
xlabel('NODE NUMBER');
ylabel('DISPLACEMENT');
 
subplot(3,2,2); 
plot(sec(:,1));
title('MODE SHAPE for 2nd natural frequency(displacement)')
xlabel('NODE NUMBER');
ylabel('DISPLACEMENT');
 
subplot(3,2,3); 
plot(third(:,1));
title('MODE SHAPE for 10th natural frequency(displacement)')
xlabel('NODE NUMBER');
ylabel('DISPLACEMENT');
 
I=zeros((n1/2)-1,1);
k=1;
for k=1:(n1/2)-1
   I(k) = I(k) + v1(2*k);
    k=k+1;
end
 
II=zeros((n1/2)-1,1);
k=1;
for k=1:(n1/2)-1
   II(k) = II(k) + v2(2*k);
    k=k+1;
end
 
III=zeros((n1/2)-1,1);
k=1;
for k=1:(n1/2)-1
    III(k) = III(k) + v3(2*k);
    k=k+1;
end
 
%MODE SHAPES for 1st three natural frequency rotational degree of freedom
subplot(3,2,4);
plot(I(:,1));
title('MODE SHAPE for 1st natural frequency(slope)')
xlabel('NODE NUMBER');
ylabel('SLOPE');
 
subplot(3,2,5); 
plot(II(:,1));
title('MODE SHAPE for 2nd natural frequency(slope)')
xlabel('NODE NUMBER');
ylabel('SLOPE');
 
subplot(3,2,6); 
plot(III(:,1));
title('MODE SHAPE for 10th natural frequency(slope)')
xlabel('NODE NUMBER');
ylabel('SLOPE');
