function plot_system(t_list, V_list, t_list_lin, V_list_lin, V_modal, t_max, tit, writerObj, box_params)
    
    % Create tiled figure
    fig1 = figure;
    fig1.Position = [0 0 1500 1000];
    t = tiledlayout(3,5);
    title(t, "Vibration Mode "+tit)

    % X Position Plots
    nexttile(1,[1,2])
    x_plots = struct();
    x_plots.lin = plot(t_list_lin, V_list_lin(:,1),".",DisplayName="Linearized", MarkerSize=3);
    hold on;
    x_plots.non = plot(0, 0,"--",DisplayName="Non-Linearized");
    x_plots.mod = plot(0, 0,"-",DisplayName="Modal");
    ylabel("X");
    xlim([0,t_max])
    legend();

    % Y Position Plots
    nexttile(6,[1,2])
    y_plots = struct();
    y_plots.lin = plot(t_list_lin, V_list_lin(:,2),".",DisplayName="Linearized", MarkerSize=3);
    hold on;
    y_plots.non = plot(0, 0,"--",DisplayName="Non-Linearized");
    y_plots.mod = plot(0, 0,"-",DisplayName="Modal");
    ylabel("Y");
    xlim([0,t_max])
    legend();

    % Theta plots
    nexttile(11,[1,2])
    t_plots = struct();
    t_plots.lin = plot(t_list_lin, V_list_lin(:,3),".",DisplayName="Linearized", MarkerSize=3);
    hold on;
    t_plots.non = plot(0, 0,"--",DisplayName="Non-Linearized");
    t_plots.mod = plot(0, 0,"-",DisplayName="Modal");
    ylabel("Theta");
    xlim([0,t_max])
    legend();
    xlabel("Time (s)")
    
    nexttile(3,[3,3])

    num_zigs = 5;
    w = .1;
    hold on;

    % Preallocate plot handles for boundary and spring drawings
    sides = [];
    for i=1:length(box_params.boundary_pts)
        sides = [sides, plot(0,0,'k','linewidth',2)];
    end
    springs = [];
    for i=1:length(box_params.k_list)
        springs = [springs, initialize_spring_plot(num_zigs,w)];
    end
    txt = text(-14,14,string(t_list(1)));
    axis equal; axis square;
    axis([-15,15,-15,15]);
    t_diff = diff(t_list); % frame timing

    % Main animation loop
    for j=1:length(t_diff)
        tic()

        % Evaluate geometry of box and environment at current state
        P1 = box_params.P_world;
        P2 = compute_rbt(V_list(j,1),V_list(j,2),V_list(j,3),box_params.P_box);
        P3 = compute_rbt(V_list(j,1),V_list(j,2),V_list(j,3),box_params.boundary_pts);
        
        % Update drawn boundary polygon
        for i=1:length(sides)
            set(sides(i),'xdata',[P3(1,i),P3(1,mod(i,length(sides))+1)],'ydata',[P3(2,i),P3(2,mod(i,length(sides))+1)]);
        end
        
        % Update each spring according to current endpoints
        for i=1:length(springs)
            update_spring_plot(springs(i),P1(:,i),P2(:,i));
        end
        
        % Update animation timestamp
        set(txt, 'string',string(t_list(j+1)))

        % Update time-history plots with newly computed values
        set(x_plots.non,'xdata',t_list(1:j+1),'ydata',V_list(1:j+1, 1));
        set(x_plots.mod,'xdata',t_list(1:j+1),'ydata',V_modal(1:j+1, 1));
        set(y_plots.non,'xdata',t_list(1:j+1),'ydata',V_list(1:j+1, 2));
        set(y_plots.mod,'xdata',t_list(1:j+1),'ydata',V_modal(1:j+1, 2));
        set(t_plots.non,'xdata',t_list(1:j+1),'ydata',V_list(1:j+1, 3));
        set(t_plots.mod,'xdata',t_list(1:j+1),'ydata',V_modal(1:j+1, 3));

        drawnow;
        %capture a frame (what is currently plotted)
        current_frame = getframe(fig1);
        %write the frame to the video
        writeVideo(writerObj,current_frame);
        p=toc();
        pause(max(0,t_diff(j)-p));
    end
end