load('vertices.dat');
load('triangles.dat');
load('boundary.dat');

nb_vertices = size(vertices,1);
nb_triangles = size(triangles,1);
nb_boundary = size(boundary,1);

% In te vullen
T = 25+273.15;
n_u = 20.8/100;
n_v = 0.04/100;


Du_r = 2.8e-10;
Du_z = 1.1e-9;
Dv_r = 2.32e-9;
Dv_z = 6.97e-9;
hu = 7e-7;
hv = 7.5e-7;
patm = 101300;
Rg = 8.314;
Cu_amb = (patm*n_u)/(Rg*T);
Cv_amb = (patm*n_v)/(Rg*T);
Vmu_ref = 2.39e-4;
E_a_vmu_ref = 80200;

Vmu = Vmu_ref*exp((E_a_vmu_ref/Rg)*((1/293.15)-(1/T)));
Kmu = 0.4103;
Kmv = 27.2438;
rq = 0.97;
Vmfv_ref = 1.61e-4;
E_a_vmfv_ref = 56700;
Vmfv = Vmfv_ref*exp((E_a_vmfv_ref/Rg)*((1/293.15)-(1/T)));
Kmfu = 0.1149;

A_u = zeros(nb_vertices,nb_vertices);
A_v = zeros(nb_vertices,nb_vertices);
B = zeros(nb_vertices,nb_vertices);
C = zeros(nb_vertices,nb_vertices);
D = zeros(nb_vertices,1);
for i=1:nb_triangles
    area = Opp_triangle(vertices(triangles(i,:),:));
    A_u(triangles(i,:),triangles(i,:)) = A_u(triangles(i,:),triangles(i,:)) + A_local(area,vertices(triangles(i,:),:),Du_r,Du_z);
    A_v(triangles(i,:),triangles(i,:)) = A_v(triangles(i,:),triangles(i,:)) + A_local(area,vertices(triangles(i,:),:),Dv_r,Dv_z);
    B(triangles(i,:),triangles(i,:)) = B(triangles(i,:),triangles(i,:)) + B_local(area,vertices(triangles(i,:),:));
end

for i=1:nb_boundary
    length = Length_edge(vertices(boundary(i,:),:));
    C(boundary(i,:),boundary(i,:)) = C(boundary(i,:),boundary(i,:))+ C_local(length,vertices(boundary(i,:),:));
    D(boundary(i,:)) = D(boundary(i,:)) + D_local(length,vertices(boundary(i,:),:));
end
F = @(u,v) [A_u*u+B*Ru(u,v,Vmu,Kmu,Kmv)+hu*(C*u-D*Cu_amb);-A_v*v+B*Rv(u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv)-hv*(C*v-D*Cv_amb)];
J = @(u,v) [[A_u+B*dRudu(u,v, Vmu,Kmu,Kmv)+hu*C B*dRudv(u,v, Vmu,Kmu,Kmv)];[B*dRvdu(u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv) -A_v+B*dRvdv(u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv)-hv*C]];

%x_0 = [(A_u+(Vmu/Kmu)*B)+hu*C zeros(nb_vertices,nb_vertices);B*rq*(Vmu/Kmu) -A_v-hv*C]\[hu*D*Cu_amb;-hv*D*Cv_amb];

%u_0 = x_0(1:nb_vertices);
%v_0 = x_0(nb_vertices+1:end);
u_0 = Cu_amb*ones(nb_vertices,1);
v_0 = Cv_amb*ones(nb_vertices,1);
[u,v] = newton_raphson( F,J,u_0,v_0,5*10^(-6));

min = -5;
max = 5;
xlin = linspace(min,max,300);
ylin = linspace(min,max,300);
[X,Y] = meshgrid(xlin,ylin);
U = griddata(vertices(:,1),vertices(:,2),u,X,Y,'linear');
V = griddata(vertices(:,1),vertices(:,2),v,X,Y,'linear');
figure
subplot(1,2,1)
contourf(X,Y,U,10)
xlim([min max])
ylim([min max])
subplot(1,2,2)
contourf(X,Y,V,10)
xlim([min max])
ylim([min max])
figure
subplot(1,2,1)
mesh(X,Y,U)
xlim([min max])
ylim([min max])
subplot(1,2,2)
mesh(X,Y,V)
xlim([min max])
ylim([min max])
