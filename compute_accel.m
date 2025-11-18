%Computes the linear and angular acceleration of the box
%given its current position and orientation
%INPUTS:
%x: current x position of the box
%y: current y position of the box
%theta: current orientation of the box
%box_params: a struct containing the parameters that describe the system
%Fields:
%box_params.m: mass of the box
%box_params.I: moment of inertia w/respect to centroid
%box_params.g: acceleration due to gravity
%box_params.k_list: list of spring stiffnesses
%box_params.l0_list: list of spring natural lengths
%box_params.P_world: 2 x n list of static mounting
% points for the spring (in the world frame) (frame)
%box_params.P_box: 2 x n list of mounting points
% for the spring (in the box frame) (corners)
%
%OUTPUTS
%ax: x acceleration of the box
%ay: y acceleration of the box
%atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)

Plist_world = compute_rbt(x,y,theta,box_params.P_box); %box corners in world frame

force_vecs = compute_spring_force(box_params.k_list,box_params.l0_list, box_params.P_world, Plist_world);
dd_r_c = -[0;box_params.g] + 1/box_params.m * sum(force_vecs,2); % accel = force/mass + gravity
ax = dd_r_c(1);
ay = dd_r_c(2);
q_i_world = Plist_world - [x;y];

q_i_world = [q_i_world; zeros(1,size(q_i_world,2))];  % offsets of spring forces about box center
force_vecs = [force_vecs; zeros(1,size(force_vecs,2))];

dd_theta = 1/box_params.I * sum(cross(q_i_world,force_vecs),2); % angular accel from torque / I
atheta = dd_theta(3);

end