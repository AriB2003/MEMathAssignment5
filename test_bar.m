% example set of system parameters for testing vibrating box simulation
% from assignment 5, page 8
LW = 10; LH = 1; LG = 3;
m = 1; Ic = (1/12)*(LHˆ2+LWˆ2);
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
box_params.l0_list = l0*ones(size(P_world,2));
box_params.P_world = P_world;
box_params.P_box = P_box;
box_params.boundary_pts = boundary_pts;