% Modificado por geem (julio 2018)
% Gráfico GMRES: norma rm en funcion del numero de ciclos

function [vec_sol, vecnormy, x]= mi_GMRES_m(A,b,m1,itermax,tol,color,print, Name_Matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crea el vector solución inicial x0
n = size(A,1);
flag = 0;
%tol=1e-06; % Parámetro de entrada
maxit = itermax; %Parámetro de entrada
restart = 0; % Cuenta los reinicios
mpv = 0; % Cuenta la generación de vectores Arnoldi
x0 = zeros(size(b,1), 1);
r0 = b - A*x0;
res(1,:)=norm(r0); % Crea un vector de resíduos, con primer elemento igual a ||r0||
logres(1,:)=(norm(r0)/res(1,1));
%iter(1,:)=restart;
m = m1;  % Parámetro de entrada. Recomendado: GMRES(30). Cambiar: m=15, m=25, m=30
% contador = 1; % No se usa, está sólo para los tests
% iter = []
% norma_de_y = zeros(1,maxit+1); % No se usa, está sólo para los tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVER
tic; % Control de tiempo de ejecucion
while flag == 0
% Aquí podría ir el controlador PID de Cuevas
%     if iter(size(iter,1),:) ~=1
%         [miter]=pdrule(m,minitial,mmin,res,iter(size(iter,1),:),mstep, mmax,alpha, delta); %cab
%         m=real(miter(1,1));
%         minitial=miter(1,2);
%     else
%         m=minitial;
%     end
    r = b - A*x0;
    beta = norm(r);
    v(:,1) = (1/beta)*r;
    mpv = mpv + 1;
    w = zeros(n,m); % Matriz auxiliar?
    h = zeros(m+1,m); % Matriz de Hessenberg
    for j = 1:m % Cómo definir m? Es la misma m del GMRES(m)? Sí
        w(:,j) = A*v(:, j);
        for i = 1:j
            h(i,j) = w(:,j)'*v(:,i);
            w(:,j) = w(:,j) - h(i,j)*v(:,i);
        end
        h(j+1,j) = norm(w(:,j));
        if h(j+1,j) == 0                % ¿¿¿¿Qué pasa si h(j+1,j) == 0????
            m = j;                      % Resetea m???
            h2 = zeros(m+1,m);          % Se crea una matriz no-cuadrada H~... 
            for k = 1:m
                h2(:,k) = h(:,k);       % Se copia la matriz, columna a columna...
            end
            h = h2;
        else
            mpv = mpv + 1;
            v(:, j+1) = w(:,j)/h(j+1,j);
        end
    end
    %H0 = h;
    g=zeros(m+1,1);
    g(1,1)=beta;

    % Transformar la matriz Hessenberg en matriz triangular superior mediante
    % "plane rotations"
    for j=1:m        
        P=eye(m+1);   
        sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2)); 
        cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2)); 
        P(j,j)=cos; 
        P(j+1,j+1)=cos;
        P(j,j+1)=sin; 
        P(j+1,j)=-sin;
        h=P*h;
        g=P*g;
    end

    R=zeros(m,m);
    G=zeros(m,1);
    V=zeros(n,m);
    for k=1:m
        G(k)=g(k);
        V(:,k)=v(:,k); % Por qué?
        for i=1:m
            R(k,i)=h(k,i);
        end
    end
    
    % y: minimizer
    minimizer=R\G;
    norm_y = norm(minimizer);
    
    % xm = x0 + V*y
    xm=x0+V*minimizer;
    
    %norma_de_y(contador)=norm(minimizer); Sólo para tests
    %contador = contador + 1; Sólo para tests

    %logres(restart+1,:)=abs(g(m+1,1));
    %iter(restart+1,:)=restart+1;
    %logres(restart+1,:)=abs(g(m+1,1))/res(1,1);
    %logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,1));
    logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,:));
    lognormy(size(logres,1)+1,:)=norm_y;
    
    % Acá la norma residual de rm es igual a la última componente de g en
    % las rotaciones de Givens (ver Proposition 6.9, eq. 6.42, Saad)
    % if (abs (g(m+1,1)))/res(1,1) <tol  || size(logres,1)==maxit    %empleando ï¿½ltima componente de g como residuo
    
    if (abs (g(m+1,1))/res(1,1)) < tol  || restart == maxit    %empleando ï¿½ltima componente de g como residuo
    %if (abs (g(s+1,1))) <tol  || restart==1000    %empleando ï¿½ltima componente de g como residuo
    %if (res(size(res,1)))/res(1,1) <tol      %empleando ï¿½ltima componente de g como residuo
         flag =1;
         %residuo = (abs (g(m+1,1)))/res(1,1);
    else
        x0=xm;                        %update and restart
        restart=restart+1;
    end
end  % END OF "WHILE FLAG"
t1 = toc; % Fin del cronómetro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if print == 1
    figure(3)
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    xlabel('Number of Restart Cycles');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
end
save('logres_GMRES_fijo.mat','logres');

vec_sol = [t1 restart m*restart];
vecnormy = lognormy;
x = xm;