function simulate_box()
    close all
    %define system parameters
    LW = 20; LH = 20; LG = 5;
    m = 0.5; Ic = (1/12)*(LH^2+LW^2);
    g = 1; k = 20; 
    k_list = [0.5*k, 0.5*k, 2*k, 5*k];
    l0 = 1.5*LG;
    boundary_pts = [7   7   8   10   10  11  14  14  15  17  17  20  21  21  21  22  22  21  20  19  18  17  13  12  11  10  10  9   8   7;
                    -13 -20 -21 -21 -24 -25 -25 -21 -20 -20 -25 -25 -24 -20 -14 -13 -10 -9  -9  -7  -6  -6  -6  -6  -7  -8  -9  -11 -11 -12];
    boundary_pts_mean = mean(boundary_pts,2);
    boundary_pts = boundary_pts-boundary_pts_mean;
    P_box = [18 18  12 12;
             -6 -6  -6 -6]-boundary_pts_mean;
    P_world = [23 18   7    12;
               -6 -1    -6  -1]-boundary_pts_mean;
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

    % Solver for equilibrium state using newton solver
    Veq = multi_newton_solver(my_rate_func2, zeros(6,1), struct());
    
    % Set initial state and time window
    x0 = 0;
    y0 = 0;
    theta0 = 0;
    vx0 = 0;
    vy0 = 0;
    vtheta0 = 0;
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    tspan = [0,6];
    
    % Jacobian for linearization
    J_approx = approximate_jacobian(my_rate_func2, Veq);
    
    % Modal analysis 
    [U_mode, omega_n] = eigs(J_approx(4:6,1:3),3);
    omega_n = sqrt(-omega_n);
    
    %define location and filename where video will be stored
    input_fname = "video.avi";
    %create a videowriter, which will write frames to the animation file
    writerObj = VideoWriter(input_fname);
    %must call open before writing any frames
    open(writerObj);
    
    % Sweep perturbation magnitudes and simulate each vibration mode
    for e=-2:0
        epsilon = 10^e;
        for w=1:3
            hz = omega_n(w,w)/(2*pi);
            Vpert = Veq + epsilon*[U_mode(:,w);0;0;0];
            
            % Run linearized simulation
            my_linear_rate = @(t_in,V_in) J_approx*(V_in-Veq);
            % Vpert = Veq + epsilon*V0;
            [t_list_lin,V_list_lin,~, ~, ~] = explicit_RK_variable_step_integration(my_linear_rate, tspan, Vpert, 0.01, rk_method("fehlberg"), 4, 10^-8);
            
            %run nonlinear simulation
            my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
            [t_list,V_list,~, ~, ~] = explicit_RK_variable_step_integration(my_rate_func, tspan, Vpert, 0.01, rk_method("fehlberg"), 4, 10^-8);
            
            % Modal analysis
            x_modal = Veq(1)+epsilon*U_mode(1,w)*cos(omega_n(w,w)*t_list);
            y_modal = Veq(2)+epsilon*U_mode(2,w)*cos(omega_n(w,w)*t_list);
            theta_modal = Veq(3)+epsilon*U_mode(3,w)*cos(omega_n(w,w)*t_list);
            
            plot_system(t_list,V_list, t_list_lin, V_list_lin, [x_modal;y_modal;theta_modal]', tspan(2), string(w)+"("+string(round(hz,2))+" Hz); epsilon="+string(round(epsilon,3)), writerObj,box_params)
        end
    end
    close(writerObj);
end