% Script para generar las matrices A y b del sistema lineal Ax = b
clear all
clc
e0 = 8.8541878176e-12; % Permitividad eléctrica del vacío
a = 1;
b = 0.25;
M = 399; 
N = 99;
k0 = 2*pi;
thita = pi/4;
Z0 = 120*pi;
hx = a/(M+1);
hy = b/(N+1);
A1= sparse(gallery('tridiag',M,1,-2,1)) .*(hx^(-2)) - eye(M).*2.*(hy^(-2));
B = (hx^(-2))*eye(M);

LIx = 0.2; LSx = 0.8; % Límite inferior y superior del subdominio en el eje x
LIy = -0.25; LSy = -0.2; % Límite inferior y superior del subdominio en el eje y
% Procedimiento para crear la matriz de permitividad
ev = zeros(N,M);
for j =1:N
    for i = 1:M
        if hx*i >= LIx && hx*i <= LSx && hy*(-j) >= LIy && hy*(-j) <= LSy
                                     % Condiciones para que (xi, yj) caiga
                                     % dentro del dominio (0.2,0.8)x(-0.25,-0.2)
            ev(j,i) = 2;
        else
            ev(j,i) = 1;
        end
    end
end

d=N; T = cell(1,d);
for j=1:d
  T{j} = sparse(A1+k0^2*diag(ev(j,:)));
end

Gr = zeros(M,M);
Gi = zeros (M,M);
t = zeros(M,M);
g = zeros(M,1);
for l=1:M
    for i=1:M
        if abs(i-l) == 1
            t(i,l) = (1-log(2))/hx; %% log: logaritmo natural
        elseif i == l
            t(i,l) = -2/hx;
        else
            t(i,l) = (log((i-l)^2/((i-l)^2 - 1)))/hx;
        end
   
        if i == l
            Gr(i,l) = -t(i,l)*k0*abs(i*hx-l*hx)*(-2/pi); % ver Jianming Jing p. 379
            Gi(i,l) = k0*k0*hx/2*(1/2); % Notar que J~(n,x)=(x/2)^n/(n!) para n>0
        else
            Gr(i,l) = -t(i,l)*k0*abs(i*hx-l*hx)*bessely(1,(k0*abs(i*hx-l*hx)))/2;
            Gi(i,l) = k0*hx*besselj(1,(k0*abs(i*hx-l*hx)))/(2*abs(i*hx-l*hx));
        end
        g(i,1) = -2*sqrt(-1)*(k0*cos(thita))*exp(sqrt(-1)*(k0*sin(thita))*i*hx);
    end    
end

G=Gr+Gi*sqrt(-1);
A2=(hy^(-1))*G-B;
lr=zeros(M,(M*(N-1)));
% Esta es la matriz A:
matrix = sparse(blkdiag(T{:})) + sparse(blktridiag(zeros(M),B,B,N));
A3 = sparse([matrix; [lr B]]);
lc=zeros(M*(N-1),M);
lc=[lc;B;A2];
A =[A3 lc];

b = zeros(M*N,1);
for count=1:M*N
    b(i)=sqrt(-1)*k0*Z0*rand;
end
b = [b;g];
% Guardar A y b en un archivo .mat    
save('cavity09.mat','A','b');