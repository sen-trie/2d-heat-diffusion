function heatDiff(test_case)
    close all; 
    clc;
    
    current_t = 0;
    max_t = 10;
    
    %checks to see if the specified corner is a neumann B.C
    %edge order: [top, bot, left, right] 
    num_placeholder = 1e-15;
    neumann_array = [num_placeholder, num_placeholder, num_placeholder, num_placeholder];
    k=10;
    [neumann_array, x, y, t_step] = setGrid(test_case, neumann_array);
    
    x_step = mean(diff(x));
    y_step = mean(diff(y));
    alpha_x = k * t_step / (x_step^2);
    alpha_y = k * t_step / (y_step^2);
    
    T=ones(length(x),length(y));
    T=boundary(test_case, T, x, y);
    T_new = T;
    figure;
    disp('Using parameters: ')
    displayBC(neumann_array, num_placeholder, T);
    disp('Running...');

    % Initialize time intervals for saving pictures (for logging purposes)
    %time_intervals = [0.1, 0.5, 10];
    %save_count = 1;
    steady_threshold = 0.01;
    steady_state = false;

    while current_t < max_t
        T_new = applyDiff(T_new, T, x, y, num_placeholder, alpha_x, alpha_y, neumann_array, x_step, y_step);
    
        % Reapply boundary conditions
        T_new=boundary(test_case, T_new, x, y);
        current_t = current_t + t_step;

        % Checking steady state
        max_temp_change = max(max(abs(T_new - T)));
        if max_temp_change <= steady_threshold
            disp(['System reached steady state at t = ', num2str(current_t), 's']);
            steady_state = true;
        end
        
        % Plotting in real-time
        pcolor(x, y, T_new), shading interp, colormap jet
        if steady_state == true
            title_text = ['System reached steady-state at time: ', num2str(current_t, '%.2f'), 's'];
        else
            title_text = ['T (^oC) at time: ', num2str(current_t, '%.2f'), 's'];
        end
        
        switch test_case
            case 'A'
                pbaspect([1 1 1])
            case 'B'
                pbaspect([1 1 1])
            case 'C'
                pbaspect([1 0.5 1])
            case 'D'
                pbaspect([1 1 1])
            case 'E'
                pbaspect([1 1 1])
            case 'F'
                pbaspect([1 1 1])
            otherwise
                error('Please insert a test case');
        end

        colorbar, xlabel('x(m)'), ylabel('y(m)'), title(title_text)
        drawnow;

        if steady_state == true
            %saveas(gcf, ['Temperature_', test_case, '_steady','.png']); % (for logging purposes)
            break;
        end

        % Save pictures at specified time intervals (for logging purposes)
        %if save_count <= length(time_intervals) && current_t >= time_intervals(save_count)
        %    saveas(gcf, ['Temperature_', test_case, '_',num2str(time_intervals(save_count)), 's','.png']);
        %    save_count = save_count + 1;
        %end

        T = T_new;
    end
end
%%
function T = boundary(tcase, T, x, y)
    switch tcase
        case 'A'
            T(end,:)=100;
            T(1,:)=100;
            T(:,1)=100;
            T(:,end)=100;
        case 'B'
            T(end,:)=40;
            T(1,:)=20;
            T(:,1)=100;
            T(:,end)=60;
        case 'C'
            T(end,:)=0;
            T(1,:)=0;
            T(:,1)=40;
            T(:,end)=80;
        case 'D'
            T(end,:)=40;
            T(1,:)=20;
            T(:,1)=100;
            T(:,end)=60;
            middle_start = round(length(x) / 2) - 1;
            middle_end = round(length(y) / 2) + 1;
            T(middle_start:middle_end, middle_start:middle_end) = 10;
        case 'E'
        case 'F'
            middle_start = round(length(x) / 2) - 2;
            middle_end = round(length(y) / 2) + 2;
            T(middle_start:middle_end, middle_start:middle_end) = 100;
        otherwise
            error('Please insert a test case');
    end
end

function [neumann_array, x, y, t_step] = setGrid(test_case, neumann_array)
    switch test_case
        case 'A'
            t_step = 0.01;
            dx = 1;
            dy = 1;
            x = 0:dx:10;
            y = 0:dy:10;
        case 'B'
            t_step = 0.001;
            dx = 0.25;
            dy = 0.25;
            x = 0:dx:10;
            y = 0:dy:10;
        case 'C'
            t_step = 0.01;
            dx = 2;
            dy = 0.5;
            x = 0:dx:20;
            y = 0:dy:5;
        case 'D'
            t_step = 0.005;
            dx = 0.5;
            dy = 0.5;
            x = 0:dx:10;
            y = 0:dy:10;
        case 'E'
            t_step = 0.01;
            dx = 1;
            dy = 1;
            x = 0:dx:10;
            y = 0:dy:10;
            neumann_array = [ -200, -700, 400, -100];
        case 'F'
            t_step = 0.01;
            dx = 1;
            dy = 1;
            x = 0:dx:10;
            y = 0:dy:10;
            neumann_array = [0, 0, 0, 0];
        otherwise
            error('Invalid test case.');
    end
end

function displayBC(neumann_array, num_placeholder, T)
    if neumann_array(4) == num_placeholder
        disp(['Top: ', num2str(T(end, end-2)), '°C (Dir B.C)']); 
    else
        disp(['Top: ', num2str(neumann_array(4)), '°C/m (Neu B.C)']); 
    end

    if neumann_array(3) == num_placeholder
        disp(['Bottom: ', num2str(T(1, 2)), '°C (Dir B.C)']); 
    else
        disp(['Bottom: ', num2str(neumann_array(3)), '°C/m (Neu B.C)']); 
    end

    if neumann_array(2) == num_placeholder
        disp(['Left: ', num2str(T(2, 1)), '°C (Dir B.C)']); 
    else
        disp(['Left: ', num2str(neumann_array(2)), '°C/m (Neu B.C)']); 
    end

    if neumann_array(1) == num_placeholder
        disp(['Right: ', num2str(T(2, end)), '°C (Dir B.C)']);
    else
        disp(['Right: ', num2str(neumann_array(1)), '°C/m (Neu B.C)']);
    end
end

function T_new = applyDiff(T_new, T, x, y, num_placeholder, alpha_x, alpha_y, neumann_array, x_step, y_step)
    for i = 1:length(x)
        for j = 1:length(y)
            if (i == 1)
                if (neumann_array(3) ~= num_placeholder)
                    if j == 1 % bottom left corner
                        T_new(i, j) = T(i, j) + alpha_x * (2*T(i+1, j) - 2*T(i, j) - 2*x_step*neumann_array(3)) + ...
                                                alpha_y * (2*T(i, j+1) - 2*T(i, j) + 2*y_step*neumann_array(2));
                    elseif j == length(y) % top left corner
                        T_new(i, j) = T(i, j) + alpha_x * (2*T(i+1, j) - 2*T(i, j) - 2*x_step*neumann_array(3)) + ...
                                                alpha_y * (2*T(i, j-1) - 2*T(i, j) + 2*y_step*neumann_array(1));
                    else
                        T_new(i, j) = T(i, j) + alpha_x * (2*T(i+1, j) - 2*T(i, j) - 2*x_step*neumann_array(3)) + ...
                                                alpha_y * (T(i, j+1) + T(i, j-1) - 2*T(i, j));
                    end
                end
            elseif (i == length(x))
                if (neumann_array(4) ~= num_placeholder)
                    if j == 1 % bottom right corner
                        T_new(i, j) = T(i, j) + alpha_x * (2*T(i-1, j) - 2*T(i, j) + 2*x_step*neumann_array(4)) + ...
                                                alpha_y * (2*T(i, j+1) - 2*T(i, j) + 2*y_step*neumann_array(2));
                    elseif j == length(y) % top right corner
                        T_new(i, j) = T(i, j) + alpha_x * (2*T(i-1, j) - 2*T(i, j) + 2*x_step*neumann_array(4)) + ...
                                                alpha_y * (2*T(i, j-1) - 2*T(i, j) + 2*y_step*neumann_array(1));
                    else 
                        T_new(i, j) = T(i, j) + alpha_x * (2*T(i-1, j) - 2*T(i, j) + 2*x_step*neumann_array(4)) + ...
                                                alpha_y * (T(i, j+1) + T(i, j-1) - 2*T(i, j));
                    end
                end
            elseif (j == 1) %Code here assumes the current point is not at x-boundary since corner cases take care of those
                if (neumann_array(2) ~= num_placeholder)
                    T_new(i, j) = T(i, j) + alpha_x * (T(i+1, j) + T(i-1, j) - 2*T(i, j)) + ...
                                            alpha_y * (2*T(i, j+1) - 2*T(i, j) - 2*y_step*neumann_array(2));
                end
            elseif (j == length(y))
                if (neumann_array(1) ~= num_placeholder)
                    T_new(i, j) = T(i, j) + alpha_x * (T(i+1, j) + T(i-1, j) - 2*T(i, j)) + ...
                                            alpha_y * (2*T(i, j-1) - 2*T(i, j) + 2*y_step*neumann_array(1));
                end
            else 
                T_new(i, j) = T(i, j) + alpha_x * (T(i+1, j) + T(i-1, j) - 2*T(i, j)) + ...
                                    alpha_y * (T(i, j+1) + T(i, j-1) - 2*T(i, j));
            end
        end
    end
end