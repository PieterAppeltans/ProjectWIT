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


Du_r = 2.8*10^(-10);
Du_z = 1.1*10^(-9);
Dv_r = 2.32*10^(-9);
Dv_z = 6.97*10^(-9);
hu = 7*10^(-7);
hv = 7.5*10^(-7);
patm = 101300;
Rg = 8.314;
Cu_amb = (patm*n_u)/(Rg*T);
Cv_amb = (patm*n_v)*(Rg*T);
Vmu_ref = 2.39*10^(-4);
E_a_vmu_ref = 80200;

Vmu = Vmu_ref*exp((E_a_vmu_ref/Rg)*(1/293.15-1/T));
Kmu = 0.4103;
Kmv = 27.2438;
rq = 0.97;
Vmfv_ref = 1.61*10^(-4);
E_a_vmfv_ref = 56700;
Vmfv = Vmfv_ref*exp((E_a_vmfv_ref/Rg)*(1/293.15-1/T));
Kmfu = 0.1149;

A_u = zeros(nb_vertices,nb_vertices);
A_v = zeros(nb_vertices,nb_vertices);
B = zeros(nb_vertices,nb_vertices);
C = zeros(nb_vertices,nb_vertices);
D = zeros(nb_vertices,1);
for i=1:nb_triangles
    area = Opp_triangle(vertices(triangles(i,:),:));
    A_u(triangles(i,:),triangles(i,:)) = A(triangles(i,:),triangles(i,:)) + A_local(area,vertices(triangles(i,:),:),Du_r,Du_z);
    A_v(triangles(i,:),triangles(i,:)) = A(triangles(i,:),triangles(i,:)) + A_local(area,vertices(triangles(i,:),:),Du_r,Du_z);
    B(triangles(i,:),triangles(i,:)) = B(triangles(i,:),triangles(i,:)) + B_local(area,vertices(triangles(i,:),:));
end

for i=1:nb_boundary
    length = Length_edge(vertices(boundary(i,:),:));
    C(boundary(i,:),boundary(i,:)) = C(boundary(i,:),boundary(i,:))+ C_local(length,vertices(boundary(i,:),:));
    D(boundary(i,:)) = D(boundary(i,:)) + D_local(length,vertices(boundary(i,:),:));
end
F = @(u,v) [A_u*u+B*Ru(u,v,Vmu,Kmu,Kmv)+hu*(C*u+D*Cu_amb);-A_v*v+B*Rv(u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv)-hv*(C*u+D*Cv_amb)];
J = @(u,v) [[A_u+B*dRudu(u,v, Vmu,Kmu,Kmv)+hu*C B*dRudv(u,v, Vmu,Kmu,Kmv)];[B*dRvdu(u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv) -A_v+B*dRvdv(u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv)-hv*C]];

%u_0 = [(Du*A+(Vmu/Kmu)*B)+hu*C]\(D*Cu_amb);
%v_0 = [(-Dv*A-hv*C)]\[-rq*B*u_0+hv*D*Cv_amb];
u_0 = 0.12*ones(nb_vertices,1);
v_0 = 0.12*ones(nb_vertices,1);
[u,v]= newton_raphson( F,J,u_0,v_0,10^(-2));
surface(vertices,u_0)