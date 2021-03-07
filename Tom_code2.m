% %Initialization 
% % placeholder = zeros(75,36);
% % delta_x = .2;
% % delta_y = sqrt(.03);
% % size_placeholder = size(placeholder);
% % x_value = 0;
% % y_value = 0;
% % nudge_index = 0;
% % counter_index = 0;
% % y_index = 0;
% % f1 = figure;
% % f2 = figure;
% % f3 = figure;
% 
% %Creates the initial nodes
% 
% for col_index = 1:size_placeholder(2)
%     
%     if y_index == 2
%         y_value = y_value + delta_y;
%         y_index = 0;
%     end
%   
%     if counter_index == 4
%         x_value = 0;
%         counter_index = 0;
%         nudge_index = 1;
%     elseif nudge_index == 2
%         x_value = delta_x / 2;
%         %nudge_index = 0;
%     end
%     
%     for row_index = 1:size_placeholder(1)
%         if rem(col_index,2) == 0 % even column so y values
%             placeholder(row_index,col_index) = y_value;
%         else
%             placeholder(row_index,col_index) = x_value;
%             x_value = delta_x + x_value;
%         end
%     end
%     
%     nudge_index = nudge_index + 1;
%     counter_index = counter_index + 1;
%     y_index = y_index + 1;
% end
% 
% %%placeholder made here
% %%%%%Start here
% 
% 
% x_col = 1;
% y_col = 2;
% for plotting = 1:18
%     hold on
%     x_section = placeholder(:,x_col);
%     y_section = placeholder(:,y_col);
%     figure(f1)
%     plot(x_section,y_section,'b.')
%     x_col = x_col + 2;
%     y_col = y_col + 2;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization 

placeholder_x = zeros(75,19);
placeholder_y = zeros(75,19);
delta_x = .2;
delta_y = sqrt(.03);
size_placeholder = size(placeholder_x);
x_value = 0;
y_value = 0;
shift_index = 0;
counter_index = 0;
y_index = 0;
f1 = figure;
f2 = figure;
f3 = figure;

x_col = 1;
y_col = 1;

for col_index = 1:size_placeholder(2)
    for row_index = 1:size_placeholder(1)
        if rem(col_index,2) == 1 % starting is 0
            placeholder_x(row_index,col_index) = x_value;
        else % starting is .1
            placeholder_x(row_index,col_index) = x_value + delta_x/2;
        end
        x_value = delta_x + x_value;
    end
    x_value = 0;
end

for col_index = 1:size_placeholder(2)
    for row_index = 1:size_placeholder(1)
        placeholder_y(row_index,col_index)= y_value;
    end
    y_value = y_value + delta_y;
end

for plotting = 1:18
    hold on
    x_section = placeholder_x(:,x_col);
    y_section = placeholder_y(:,y_col);
    figure(1)
    plot(x_section,y_section,'b.')
    x_col = x_col + 1;
    y_col = y_col + 1;
end

%
%insert code for cuffs
%

cuffWidthHeight = 2.95;
cuffStart = 1.5;
for cuffWidth = 1:13
    xCuff(cuffWidth) = cuffStart;
    cuffStart = cuffStart + 1;
    yCuff(cuffWidth) = 2.95;
    plot(xCuff, yCuff, 'k-', 'LineWidth', 1.5);
end
hold off
%initialization for compression and nudging

compress_x = zeros(75,19);
compress_y = zeros(75,19);
d_o = delta_x;
index_n = 0; %x-position on the graph
index_m = 0; %y-position on the graph
x_val = 0; %x-position value; temporary storage 
y_val = 0; %y-position value; temporary storage


%placeholder will be input for this function
%Iterate x and y positions from placeholder



% Net_Fx > 0 -> shitft = +0.0005
%If the net x-directed force is > 0 , then nudge the node 5 microns to the right, etc.  If the net x-directed force is < 0 , then nudge the node 5 microns to the left, and similarly for the y-dimension.  
% for possible = 1:500. 
    %per 1 iteration, do we have to test each node's net force?(75*36)
    % per each iteration, we find the positions of all nodes and then
    % calculate the netfore for each node and decide whether they should
    % shift right or left? 
    
%force equation
%Initialization
% x_ne = 0;
% x_nw = 0 ;
% x_e = 0;
% x_w = 0;
% x_sw = 0;
% x_se = 0;



%To index from placeholder[inputarray]
%indexing x-position(x,2y-1)
%indexing y-position(x,2y)

%Compression
%iterates through all the nodes



for index_n = 1:75 %x
    for index_m = 1:19  %y
        %find main node; the node of interest; we will find nw, ne, sw, se
        %based on this node
        %Find node of interest
        x_val = placeholder_x(index_n, index_m);
        y_val = placeholder_y(index_n, index_m);
        if abs(7.5 - x_val) <= 6 
            compress_x(index_n, index_m) = x_val;
            compress_y(index_n, index_m) = y_val * 0.9;
        else
            compress_x(index_n, index_m) = x_val;
            compress_y(index_n, index_m) = y_val;
        end
    end
end

x_col = 1;
y_col = 1;
for plotting_23 = 1:18
    hold on
    x_section = compress_x(:,x_col);
    y_section = compress_y(:,y_col);
    figure(f2)
    plot(x_section,y_section,'b.')
    x_col = x_col + 1;
    y_col = y_col + 1;
end
hold off
%k

nudge_x = compress_x;
nudge_y = compress_y;
nudge_temp_x = zeros(75,19);
nudge_temp_y = zeros(75,19);
k = 1; %k eventualy cancels out during net force calculations. 
index_force = 0;
while index_force < 500
    for index_n = 1:75 %y   m - 1 1:18
        for index_m = 1:19 %x
            %find main node; the node of interest; we will find nw, ne, sw, se
            %based on this node
            %Find node of interest
            x_val = nudge_x(index_n, index_m);
            y_val = nudge_y(index_n, index_m);

            %use index_m to determine if it is odd row or even row
            %for instance, m represent 1,2,3,4,5,6,7,8,9th row and so on.
            %depending on these row's odd or evenness, the positions of nw, ne
            %and so on

            %find adjacent node positions
            %we will always code in this order
            %nw -> ne -> w-> e -> sw -> se
            if rem(index_m, 2) == 1 %odd index_m but with a EVEN ROWS
                %nw
                %index_n -1 , index_m + 1
                
                
                if index_m == 19 || index_n == 1
                    x_nw = x_val; %making the equation zero
                    y_nw = y_val;
                else
                    x_nw = nudge_x(index_n - 1, index_m + 1);
                    y_nw = nudge_y(index_n - 1, index_m + 1);
                end
                
                %ne
                %index_n, index_m +1
                
                if index_m == 19 || index_n == 75
                    x_ne = x_val; %making the equation zero
                    y_ne = y_val;
                else
                    x_ne = nudge_x(index_n, index_m + 1);
                    y_ne = nudge_y(index_n, index_m + 1);
                end

                %w
                %index_n - 1, index_m
                
                if index_n == 1
                    x_w = x_val; %making the equation zero
                    y_w = y_val;
                else
                    x_w = nudge_x(index_n - 1, index_m);
                    y_w = nudge_y(index_n - 1, index_m);
                end
                
                %e
                %index_n + 1, index_m
                
                if index_n == 75
                    x_e = x_val; %making the equation zero
                    y_e = y_val;
                else
                    x_e = nudge_x(index_n + 1, index_m);
                    y_e = nudge_y(index_n + 1, index_m); 
                end
                
                %sw 
                %index_n - 1 , index_m - 1
                
                if index_m == 1 || index_n == 1
                    x_sw = x_val; %making the equation zero
                    y_sw = y_val;
                else
                    x_sw = nudge_x(index_n - 1, index_m - 1);
                    y_sw = nudge_y(index_n - 1, index_m - 1);

                end
                
                %se 
                %index_n, index_m - 1
                
                if index_m == 1 || index_n == 75
                    x_se = x_val; %making the equation zero
                    y_se = y_val;
                else
                    x_se = nudge_x(index_n, index_m - 1);
                    y_se = nudge_y(index_n, index_m - 1);
                end
                
            else    %even *******************************************
                %nw
                %index_n, index_m + 1
                
                
                if index_m == 19 || index_n == 1
                    x_nw = x_val; %making the equation zero
                    y_nw = y_val;
                else
                    x_nw = nudge_x(index_n , index_m + 1);
                    y_nw = nudge_y(index_n , index_m + 1);
                end
                
                %ne
                %index_n + 1, index_m +1
                
                
                if index_m == 19 || index_n == 75
                    x_ne = x_val; %making the equation zero
                    y_ne = y_val;
                else
                    x_ne = nudge_x(index_n + 1, index_m + 1);
                    y_ne = nudge_y(index_n + 1, index_m + 1);
                end
                
                %w
                %index_n - 1, index_m
                
                if  index_n == 1
                    x_w = x_val; %making the equation zero
                    y_w = y_val;
                else
                    x_w = nudge_x(index_n - 1, index_m);
                    y_w = nudge_y(index_n - 1, index_m);
                end

                %e
                %index_n + 1, index_m
                
                if  index_n == 75
                    x_e = x_val; %making the equation zero
                    y_e = y_val;
                else
                    x_e = nudge_x(index_n + 1, index_m);
                    y_e = nudge_y(index_n + 1, index_m);
                end
                
                %sw 
                %index_n, index_m - 1
                 
                if index_m == 1 || index_n == 1
                    x_sw = x_val; %making the equation zero
                    y_sw = y_val;
                else
                    x_sw = nudge_x(index_n, index_m - 1);
                    y_sw = nudge_y(index_n, index_m - 1);
                end
                
                %se 
                %index_n + 1, index_m - 1
                
                if index_m == 1 || index_n == 75
                    x_se = x_val; %making the equation zero
                    y_se = y_val;
                else
                    x_se = nudge_x(index_n + 1, index_m - 1);
                    y_se = nudge_y(index_n + 1, index_m - 1);
                end
            end

            
            
            %x
            %%Force Equations enter
            dnw = sqrt((x_nw - x_val)^2 + (y_nw - y_val)^2);
            dne = sqrt((x_ne - x_val)^2 + (y_ne - y_val)^2);
            dw = sqrt((x_w - x_val)^2 + (y_w - y_val)^2);
            de = sqrt((x_e - x_val)^2 + (y_e - y_val)^2);
            dsw = sqrt((x_sw - x_val)^2 + (y_sw - y_val)^2);
            dse = sqrt((x_se - x_val)^2 + (y_se - y_val)^2); 

            fx_ne = k * (dne - delta_x) * ((x_ne - x_val)/dne);
            fx_nw = k * (dnw - delta_x) * ((x_nw - x_val)/dnw);
            fx_se = k * (dse - delta_x) * ((x_se - x_val)/dse);
            fx_sw = k * (dsw - delta_x) * ((x_sw - x_val)/dsw);
            fx_e = k * (de - delta_x) * ((x_e - x_val)/de);
            fx_w = k * (dw - delta_x) * ((x_w - x_val)/dw);
            
            fy_ne = k * (dne - delta_x) * ((y_ne - y_val)/dne);
            fy_nw = k * (dnw - delta_x) * ((y_nw - y_val)/dnw);
            fy_se = k * (dse - delta_x) * ((y_se - y_val)/dse);
            fy_sw = k * (dsw - delta_x) * ((y_sw - y_val)/dsw);
            fy_e = k * (de - delta_x) * ((y_e - y_val)/de);
            fy_w = k * (dw - delta_x) * ((y_w - y_val)/dw);
            
            
            net_force_x = fx_ne + fx_nw + fx_se + fx_sw + fx_e + fx_w;
            net_force_y = fy_ne + fy_nw + fy_se + fy_sw + fy_e + fy_w;
%             %directional k-prime values
%             kp_ne = k*(dne - d_o) / dne;
%             kp_nw = k*(dnw - d_o) / dnw;
%             kp_se = k*(dse - d_o) / dse;
%             kp_sw = k*(dsw - d_o) / dsw;
%             kp_e = k*(de - d_o) / de;
%             kp_w = k*(dw - d_o) / dw;

%NEED TO UPDATE FORCE EQUATION WITH THE ABOVE VARIABLES
%                 net_force_x = (kp_ne * (x_ne - x_val)) - (kp_nw* (x_nw - x_val)) + (kp_se * (x_se - x_val)) - (kp_sw * (x_sw - x_val)) +  (kp_e * (x_e - x_val)) - (kp_w * (x_w - x_val));
%                 net_force_y = (kp_ne * (y_ne - y_val)) - (kp_nw* (y_nw - y_val)) + (kp_se * (y_se - y_val)) - (kp_sw * (y_sw - y_val)) +  (kp_e * (y_e - y_val)) - (kp_w * (y_w - y_val));
            
%             net_force_x = -(kp_ne * (x_ne - x_val)) + (kp_nw* (x_nw - x_val)) - (kp_se * (x_se - x_val)) + (kp_sw * (x_sw - x_val)) - (kp_e * (x_e - x_val)) + (kp_w * (x_w - x_val));
%             net_force_y = -(kp_ne * (y_ne - y_val)) + (kp_nw* (y_nw - y_val)) - (kp_se * (y_se - y_val)) + (kp_sw * (y_sw - y_val)) -  (kp_e * (y_e - y_val)) + (kp_w * (y_w - y_val));
%             
            if net_force_x > 0 
                nudge_temp_x(index_n, index_m) = x_val + 0.0005;

            elseif net_force_x < 0 
                nudge_temp_x(index_n, index_m) = x_val - 0.0005;

            else
                nudge_temp_x(index_n, index_m) = x_val;

            end

            %y
            if net_force_y > 0 
                nudge_temp_y(index_n, index_m) = y_val + 0.0005;

            elseif net_force_y < 0 
                nudge_temp_y(index_n, index_m) = y_val - 0.0005;

            else
                nudge_temp_y(index_n, index_m) = y_val;

            end
            
            
        end
    end
   
    nudge_x = nudge_temp_x;
    nudge_y = nudge_temp_y;
    %change index number and add plotting 
    index_force = index_force + 1;
end


x_col = 1;
y_col = 1;
for plotting_2 = 1:18
    hold on
    x_section = nudge_x(:,x_col);
    y_section = nudge_y(:,y_col);
    figure(f3)
    plot(x_section,y_section,'b.-')
    x_col = x_col + 1;
    y_col = y_col + 1;
end
%k
hold off