clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'cavity07.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Name_Matrix = 'cavity07';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only for UFlorida Matrix Market
% A = Problem.A;
% b = Problem.b(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parámetros globales
itermax = 500; %2000 funciona
tol = 1e-06;
end_count = 1;
% Parámetro de Adaptive-LGMRESE     
eps0 = 0.5;
alpha = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LGMRES-E(m,d,l)
m0 = 30;
d0 = 2;
l0 = 2;
color = 'r-';
print = 0;
temp = zeros(end_count, 1);
for count = 1:end_count
    if count == end_count
        print = 1;
    end
    [vec5, vny, xx] = mi_LGMRESE(A, b, m0, d0, l0, itermax, tol, alpha, eps0, color, print, Name_Matrix);
    temp(count)=vec5(1);
end
sum_s= vec5(3);
t_prom5 = mean(temp);
t_std_dev5 = std(temp);
restart = vec5(2);
metrics_LGMRESE = [t_prom5 t_std_dev5 vec5(2) vec5(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GMRES(m)
color='k-';
m = 30; % El más recomendado: GMRES(30). Cambiar: m=15, m=25, m=30
print=0;
temp=zeros(end_count,1);
%fig_GMRES = figure(1);
for count = 1:end_count
   if count==end_count
       print = 1;
   end
   [vec1, vny, xx]= mi_GMRES_m(A,b,m,itermax,tol, color,print, Name_Matrix);
                                    % Agregar la tolerancia tol=1e-10
   temp(count)=vec1(1);
end
t_prom1=mean(temp);
t_std_dev1=std(temp);
metrics_GMRES = [t_prom1 t_std_dev1 vec1(2) vec1(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GMRES(m_k), m variable
color='b-';
m1 = 30; % El más recomendado: GMRES(30). Cambiar: m=15, m=25, m=30
print=0;
temp=zeros(end_count,1);
for count = 1:end_count
   if count==end_count
       print = 1;
   end
   [vec6]= mi_GMRES_mj(A,b,m1,itermax,tol, eps0, alpha, color,Name_Matrix,print);
                                    % Agregar la tolerancia tol=1e-10
   temp(count)=vec6(1);
end
t_prom6=mean(temp);
t_std_dev6=std(temp);
metrics_GMRES_mj = [t_prom6 t_std_dev6 vec6(2) vec6(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LGMRES(m,l) 
color='g-';
m1=27;
l1=3;
temp=zeros(end_count,1);
print=0;
for count = 1:end_count
   if count==end_count
       print = 1;
   end
   [vec2]= mi_LGMRES(A, b, m1, l1, itermax, tol, color, print, Name_Matrix);
                                    % Agregar la tolerancia tol=1e-10
   temp(count)=vec2(1,1);
end

t_prom2=mean(temp);
t_std_dev2=std(temp);
metrics_LGMRES = [t_prom2 t_std_dev2 vec2(2) vec2(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GMRES-E(m,d)
color='m-';
print=0;
mE=27;
dE=3;
temp=zeros(end_count,1);
for count = 1:end_count
   if count==end_count
       print = 1;
   end
   [vec4]= mi_GMRES_E(A, b, mE, dE, itermax, tol, color, print, Name_Matrix);
                                    % Agregar la tolerancia tol=1e-10
   temp(count)=vec4(1,1);
end
t_prom4=mean(temp);
t_std_dev4=std(temp);
metrics_GMRESE = [t_prom4 t_std_dev4 vec4(2) vec4(3)]; %revisar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend(['A-LGMRES-E(mj,',num2str(d0),',',num2str(l0),'), t=',num2str(t_prom5),'±', num2str(t_std_dev5), 's'],...
     ['GMRES(',num2str(m),'),t=', num2str(t_prom1),'±', num2str(t_std_dev1), 's'],...
     ['GMRES(mj),t=', num2str(t_prom6),'±', num2str(t_std_dev6), 's'],...
     ['LGMRES(',num2str(m1),',',num2str(l1),'), t=',num2str(t_prom2),'±', num2str(t_std_dev2), 's'],...
     ['GMRES-E(',num2str(mE),',',num2str(dE),'), t=',num2str(t_prom4),'±', num2str(t_std_dev4), 's'],...
     'Location','Northeast');

my_plots_path = 'C:\Users\HP\Documents\GitHub\mi-gmres-results';
% legend(['GMRES(',num2str(m),'),t=', num2str(t_prom1),'±', num2str(t_std_dev1), 's'],...
%      'Location','Northeast');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
% save('metrics_wang4_exp01_20181024.mat','metrics_LGMRESE','metrics_GMRES','metrics_GMRES_mj',...
%     'metrics_LGMRES','metrics_GMRESE','eps0', 'alpha', 'end_count', 'itermax','Name_Matrix');
% 
% save('metrics_sherman3_20180930.mat','metrics_GMRES', 'end_count');