function plot_system(t_list, V_list, box_params)
    figure;
    num_zigs = 5;
    w = .1;
    hold on;
    sides = [];
    for i=1:length(box_params.boundary_pts)
        sides = [sides, plot(0,0,'k','linewidth',2)];
    end
    springs = [];
    for i=1:length(box_params.k_list)
        springs = [springs, initialize_spring_plot(num_zigs,w)];
    end
    txt = text(-9,9,string(t_list(1)));
    axis equal; axis square;
    axis([-10,10,-10,10]);
    t_diff = diff(t_list);
    for j=1:length(t_diff)
        tic()
        P1 = box_params.P_world;
        P2 = compute_rbt(V_list(j,1),V_list(j,2),V_list(j,3),box_params.P_box);
        P3 = compute_rbt(V_list(j,1),V_list(j,2),V_list(j,3),box_params.boundary_pts);
        for i=1:length(sides)
            set(sides(i),'xdata',[P3(1,i),P3(1,mod(i,length(sides))+1)],'ydata',[P3(2,i),P3(2,mod(i,length(sides))+1)]);
        end
        for i=1:length(springs)
            update_spring_plot(springs(i),P1(:,i),P2(:,i));
        end
        set(txt, 'string',string(t_list(j+1)))
        drawnow;
        p=toc();
        pause(max(0,t_diff(j)-p));
    end
end