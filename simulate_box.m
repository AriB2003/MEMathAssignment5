function simulate_box()
close all
%define system parameters
LW = 10; LH = 1; LG = 3;
m = 1; Ic = (1/12)*(LH^2+LW^2);
g = 1; k = 20; k_list = [.5*k,.5*k,2*k,5*k];
l0 = 1.5*LG;
Pbl_box = [-LW;-LH]/2;
Pbr_box = [LW;-LH]/2;
Ptl_box = [-LW;LH]/2;
Ptr_box = [LW;LH]/2;
boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];
Pbl1_world = Pbl_box + [-LG;-LG];
Pbl2_world = Pbl_box + [LG;-LG];
Pbr1_world = Pbr_box + [0;-l0];
Pbr2_world = Pbr_box + [l0;0];
P_world = [Pbl1_world,Pbl2_world,Pbr1_world,Pbr2_world];
P_box = [Pbl_box,Pbl_box,Pbr_box,Pbr_box];
%define system parameters
box_params = struct();
box_params.m = m;
box_params.I = Ic;
box_params.g = g;
box_params.k_list = k_list;
box_params.l0_list = l0*ones(1,size(P_world,2));
box_params.P_world = P_world;
box_params.P_box = P_box;
box_params.boundary_pts = boundary_pts;
%load the system parameters into the rate function
%via an anonymous function
my_rate_func2 = @(V_in) box_rate_func(0,V_in,box_params);
Veq = multi_newton_solver(my_rate_func2, zeros(6,1), struct());

x0 = 1;
y0 = 1;
theta0 = 1;
vx0 = 0;
vy0 = 0;
vtheta0 = 0;
V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
tspan = [0,10];

epsilon = 10^-1;
J_approx = approximate_jacobian(my_rate_func2, Veq);

[U_mode, omega_n] = eigs(J_approx(4:6,1:3),3);
omega_n = sqrt(-omega_n);
hz = omega_n(3,3)/(2*pi)
Vpert = Veq + epsilon*[U_mode(:,3);0;0;0];

my_linear_rate = @(t_in,V_in) J_approx*(V_in-Veq);
% Vpert = Veq + epsilon*V0;
[t_list_lin,V_list_lin,~, ~, ~] = explicit_RK_variable_step_integration(my_linear_rate, tspan, Vpert, 0.01, rk_method("dormandprince"), 4, 10^-8);

%run the integration
my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
[t_list,V_list,~, ~, ~] = explicit_RK_variable_step_integration(my_rate_func, tspan, Vpert, 0.01, rk_method("dormandprince"), 4, 10^-8);

subplot(3,1,1);
plot(t_list_lin, V_list_lin(:,1),".",DisplayName="Linearized", MarkerSize=3)
hold on
subplot(3,1,2);
plot(t_list_lin, V_list_lin(:,2),".",DisplayName="Linearized", MarkerSize=3)
hold on
subplot(3,1,3);
plot(t_list_lin, V_list_lin(:,3),".",DisplayName="Linearized", MarkerSize=3)
hold on

subplot(3,1,1);
plot(t_list, V_list(:,1),"--",DisplayName="Non-Linearized")
subplot(3,1,2);
plot(t_list, V_list(:,2),"--",DisplayName="Non-Linearized")
subplot(3,1,3);
plot(t_list, V_list(:,3),"--",DisplayName="Non-Linearized")

subplot(3,1,1);
x_modal = Veq(1)+epsilon*U_mode(1,3)*cos(omega_n(3,3)*t_list);
plot(t_list, x_modal,"-",DisplayName="Modal")
ylabel("X");
legend();
subplot(3,1,2);
y_modal = Veq(2)+epsilon*U_mode(2,3)*cos(omega_n(3,3)*t_list);
plot(t_list, y_modal,"-",DisplayName="Modal")
ylabel("Y");
legend();
subplot(3,1,3);
theta_modal = Veq(3)+epsilon*U_mode(3,3)*cos(omega_n(3,3)*t_list);
plot(t_list, theta_modal,"-",DisplayName="Modal")
xlabel("Time (s)")
ylabel("Theta");
legend();

plot_system(t_list,V_list,box_params)
end