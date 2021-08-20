% Modificado por geem (febrero 2018)

function [vec_sol]= mi_GMRES_mj(A,b,m1,itermax,tol, eps0, alpha, color,Name_Matrix,print)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crea el vector solución inicial x0
n = size(A,1);
flag=0;
%tol=1e-10; % Parámtro de entrada
maxit = itermax; % Parámetro de entrada
% norma_de_y = zeros(maxit+1); % OJO!
restart=0; % Cuenta los reinicios
mpv=0; % Cuenta la generación de vectores Arnoldi
x0 = zeros(size(b,1), 1);
r0 = b - A*x0;
res(1,:)=norm(r0);
logres(1,:)=(norm(r0)/res(1,1));
%iter(1,:)=restart;
m = m1;  % El más recomendado: GMRES(30). Cambiar: m=15, m=25, m=30
m_max=100;
% contador = 1;
% iter = []
miteracion(restart+1,1)=m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLADOR
%eps0 = 1e-10; % Parámetro del controlador
tic; % Control de tiempo de ejecucion
while flag == 0
    % VARIAR m en función del minimizer "y"
    % Se busca una estrategia de control
    %           m_j = m_(j-1) + u_j
    % Donde:
    % 1) Si |y| < eps
    %       u_j = alpha
    % 2) Si no
    %       u_j = 0
    if restart > 0 && norma_de_y(restart,1) <= eps0
        if m < m_max
            m = m + alpha;
        else
            m = m_max;
        end
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% SOLVER    
    r = b - A*x0;
    beta = norm(r);
    v(:,1) = (1/beta)*r;
    mpv = mpv + 1;
    w = zeros(n,m); % OJO!
    h = zeros(m+1,m+1); % OJO!
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

    minimizer=R\G;

    xm=x0+V*minimizer;
    
    norma_de_y(restart+1,1)=norm(minimizer); %Sólo para tests

    logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,:));
    
    % Acá la norma residual de rm es igual a la última componente de g en
    % las rotaciones de Givens (ver Proposition 6.9, eq. 6.42, Saad)
    % if (abs (g(m+1,1)))/res(1,1) <tol  || size(logres,1)==maxit    %empleando ï¿½ltima componente de g como residuo
    
    if (abs (g(m+1,1))/res(1,1)) < tol  || restart==maxit    %empleando ï¿½ltima componente de g como residuo
         flag=1;
         %residuo= (abs (g(m+1,1)))/res(1,1);
    else
        x0=xm;                        %update and restart
        miteracion(restart+1,1)=m;
        restart=restart+1;
    end
end  % END OF "WHILE FLAG"
t1 = toc; % Fin del cronómetro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if print == 1
    figure(1)
    subplot(2,1,1);
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    xlabel('Number of Restart Cycles');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    subplot(2,1,2);
    semilogy(norma_de_y,color)
    ylabel('||y_j||');
    xlabel('Number of Restart Cycles')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    figure(2)
    subplot(2,1,1);
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    xlabel('Number of Restart Cycles');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    subplot(2,1,2);
    plot(miteracion, color)
    ylabel('m_j (s_j)');
    xlabel('Number of Restart Cycles')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    figure(3)
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    xlabel('Number of Restart Cycles');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    figure(4)
    subplot(2,1,1);
    plot(miteracion, color)
    ylabel('m_j (s_j)');
    xlabel('Number of Restart Cycles')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    subplot(2,1,2);
    semilogy(norma_de_y,color)
    ylabel('||y_j||');
    xlabel('Number of Restart Cycles')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on 
    
    figure(5)
    subplot(3,1,1);
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    %xlabel('Number of Restart Cycles');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    subplot(3,1,2);
    plot(miteracion, color)
    ylabel('m_j (s_j)');
    %xlabel('Number of Restart Cycles')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    subplot(3,1,3)
    semilogy(norma_de_y,color)
    ylabel('||y_j||');
    xlabel('Number of Restart Cycles')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
end
sum_m = sum(miteracion);
save('logres_GMRES_variable.mat','logres');
vec_sol = [t1 restart sum_m];