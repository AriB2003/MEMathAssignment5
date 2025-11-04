function plot_system(t_list, V_list, box_params)
    figure;
    num_zigs = 5;
    w = .1;
    hold on;
    springs = [];
    for i=1:length(box_params.k_list)
        springs = [springs, initialize_spring_plot(num_zigs,w)];
    end
    axis equal; axis square;
    axis([-10,10,-10,10]);
    t_diff = diff(t_list);
    for i=1:length(t_diff)
        P1 = box_params.P_world;
        P2 = compute_rbt(V_list(i,1),V_list(i,2),V_list(i,3),box_params.P_box);
        for i=1:length(springs)
            update_spring_plot(springs(i),P1(:,i),P2(:,i));
        end
        drawnow;
        pause(t_diff(i));
    end
end