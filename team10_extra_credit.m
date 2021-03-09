%Program Description
%BME 256 Accuracy of blood pressure readings
%Team 10
%Tom Lee
%Denny Guo
%Chloe Matasovsky
%Shwath Kumaravel
%This program demonstrates finite element modeling of the cuffed arm model
%based on cuff width and cuff indentation on the brachial artery

%Initialization 

placeholder_x = zeros(75,19); %Empty Matrix
placeholder_y = zeros(75,19); %Empty Matrix
delta_x = .2; %horizontal distance between initial nodes(cm)
delta_y = sqrt(.03); %Vertical distance between initial nodes(cm)
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


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% The following three for loops plot initial positions of the nodes
% placeholder_x -> contains x values (75x18)
% placeholder_y -> contains y values (75x18)

for col_index = 1:size_placeholder(2)
    for row_index = 1:size_placeholder(1)
        if rem(col_index,2) == 1 % starting is 0
            placeholder_x(row_index,col_index) = x_value;
        else % starting is .1cm
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
    plot(x_section,y_section,'b.-')
    x_col = x_col + 1;
    y_col = y_col + 1;
end
title("Initial Positions of the cuff and arm model(cm)");
xlabel("Length(cm)");
ylabel("Height(cm)")

%Initialization 
%Parameters for cuff size and position

model_width = 14.9;
cuffWidthHeight = 2.95;
model_midpoint = 15 / 2;
cuffwidth = 12;
cuff_diff = cuffwidth / 2;

cuffstart = model_midpoint - cuff_diff;
cuffend = model_midpoint + cuff_diff;


%Plots the cuff
cuff_array = [cuffstart, cuffend];
cuff_height = [cuffWidthHeight, cuffWidthHeight];
plot(cuff_array, cuff_height, 'k-', 'LineWidth', 1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
%Compression parameters are initialized here

compress_x = zeros(75,19);
compress_y = zeros(75,19);
d_o = delta_x;
index_n = 0; %x-position on the graph(cm)
index_m = 0; %y-position on the graph (cm)
x_val = 0; %x-position value; temporary storage; Different from x_value
y_val = 0; %y-position value; temporary storage
indent = 0.3; %Length of indentation (cm)
model_height = 3; %Height of the arm model(cm)
compressionfactor = 1 - (indent / model_height); %compression factor for the model


%This for loop implements the indentation provided into the initial model
%The calculation is done by usign absolute value between the difference
%between the model width and cuff width

for index_n = 1:75 %x
    for index_m = 1:19 %y
        %find main node; the node of interest; we will find nw, ne, sw, se
        %based on this node
        %Find node of interest
        x_val = placeholder_x(index_n, index_m);
        y_val = placeholder_y(index_n, index_m);
        if abs((model_width / 2) - x_val) <= ((cuffwidth) / 2)   %15 -> model width , 12 -> cuff width
            compress_x(index_n, index_m) = x_val;
            compress_y(index_n, index_m) = y_val * compressionfactor;
        else
            compress_x(index_n, index_m) = x_val;
            compress_y(index_n, index_m) = y_val;
        end
    end
end

x_col = 1;
y_col = 1; 
%This for loop plots the compressed nodes

for plotting_23 = 1:18
    hold on
    x_section = compress_x(:,x_col);
    y_section = compress_y(:,y_col);
    figure(f2)
    plot(x_section,y_section,'b.-')
    x_col = x_col + 1;
    y_col = y_col + 1;
end
title({['Compressed Positions of the cuff and arm model(cm)'],['[Pre-nudge: ',num2str(indent),'cm]']});
xlabel("Length(cm)");
ylabel("Height(cm)")
hold off


%%%%%%%%%%%%Initialization

row_top = 10; %Top row of the test artery
row_bot = 9; %Bottom side of the test artery
nudgenode_x = compress_x; %Transfer matrix of x values
nudgenode_y = compress_y; %Transfer matrix of y values
nudge_temp_x = zeros(75,19);
nudge_temp_y = zeros(75,19);
 %Initial Spring constant in quasiequilibrium state
%k constant eventualy cancels out during net force calculations. 
index_force = 0; %Controls the number of updates of nudging according to net force
%Adjacent nodes initialized
xadj_ne = 0;
xadj_nw = 0;
xadj_w = 0;
xadj_e = 0;
xadj_se = 0;
xadj_sw = 0;
yadj_ne = 0;
yadj_nw = 0;
yadj_w = 0;
yadj_e = 0;
yadj_se = 0;
yadj_sw = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%NUDGE NODES BASED ON NET FORCE

while index_force < 500
    for index_n = 1:75 %y   m - 1 1:18
        for index_m = 1:19 %x
            %find main node; the node of interest; we will find nw, ne, sw, se
            %based on this node
            %Find node of interest
            if index_m <= 13
                k = 10;
            else
                k = 1;
            end
            x_val = nudgenode_x(index_n, index_m);
            y_val = nudgenode_y(index_n, index_m);

            if index_m == 12
                disp(k)
            end
            %use index_m to determine if it is odd row or even row
            %for instance, m represent 1,2,3,4,5,6,7,8,9th row and so on.
            %depending on these row's odd or evenness, the positions of nw, ne
            %and so on

            %find adjacent node positions
            %we will always code in this order
            %nw -> ne -> w-> e -> sw -> se
            if rem(index_m, 2) == 1 %odd
                %nw
                %index_n -1 , index_m + 1
               % (index_m == 18 && abs((model_width / 2) - nudge_x(index_m, index_n)) > (model_width - cuffwidth))
                
                if index_n == 1 || index_m == 19
                    xadj_nw = x_val; %making the equation zero
                    yadj_nw = y_val;
                else
                    xadj_nw = nudgenode_x(index_n - 1, index_m + 1);
                    yadj_nw = nudgenode_y(index_n - 1, index_m + 1);
                    
                end
                
                %ne
                %index_n, index_m +1
                
                if index_m == 19 || index_n == 75
                    xadj_ne = x_val; %making the equation zero
                    yadj_ne = y_val;
                else
                    xadj_ne = nudgenode_x(index_n, index_m + 1);
                    yadj_ne = nudgenode_y(index_n, index_m + 1);
                end

                %w
                %index_n - 1, index_m
                
                if index_n == 1
                    xadj_w = x_val; %making the equation zero
                    yadj_w = y_val;
                else
                    xadj_w = nudgenode_x(index_n - 1, index_m);
                    yadj_w = nudgenode_y(index_n - 1, index_m);
                end
                
                %e
                %index_n + 1, index_m
                
                if index_n == 75
                    xadj_e = x_val; %making the equation zero
                    yadj_e = y_val;
                else
                    xadj_e = nudgenode_x(index_n + 1, index_m);
                    yadj_e = nudgenode_y(index_n + 1, index_m); 
                end
                
                %sw 
                %index_n - 1 , index_m - 1
                
                if index_m == 1 || index_n == 1
                    xadj_sw = x_val; %making the equation zero
                    yadj_sw = y_val;
                else
                    xadj_sw = nudgenode_x(index_n - 1, index_m - 1);
                    yadj_sw = nudgenode_y(index_n - 1, index_m - 1);
                end
                
                %se 
                %index_n, index_m - 1
                
                if index_m == 1 || index_n == 75
                    xadj_se = x_val; %making the equation zero
                    yadj_se = y_val;
                else
                    xadj_se = nudgenode_x(index_n, index_m - 1);
                    yadj_se = nudgenode_y(index_n, index_m - 1);
                end
                
            else    %even
                %nw
                %index_n, index_m + 1
                
                
                if index_m == 19 || index_n == 1
                    xadj_nw = x_val; %making the equation zero
                    yadj_nw = y_val;
                else
                    xadj_nw = nudgenode_x(index_n , index_m + 1);
                    yadj_nw = nudgenode_y(index_n , index_m + 1);
                end
                
                %ne
                %index_n + 1, index_m +1
                
                
                if index_m == 19 || index_n == 75
                    xadj_ne = x_val; %making the equation zero
                    yadj_ne = y_val;
                else
                    xadj_ne = nudgenode_x(index_n + 1, index_m + 1);
                    yadj_ne = nudgenode_y(index_n + 1, index_m + 1);
                end
                
                %w
                %index_n - 1, index_m
                
                if  index_n == 1
                    xadj_w = x_val; %making the equation zero
                    yadj_w = y_val;
                else
                    xadj_w = nudgenode_x(index_n - 1, index_m);
                    yadj_w = nudgenode_y(index_n - 1, index_m);
                end

                %e
                %index_n + 1, index_m
                
                if  index_n == 75
                    xadj_e = x_val; %making the equation zero
                    yadj_e = y_val;
                else
                    xadj_e = nudgenode_x(index_n + 1, index_m);
                    yadj_e = nudgenode_y(index_n + 1, index_m);
                end
                
                %sw 
                %index_n, index_m - 1
                 
                if index_m == 1 || index_n == 1
                    xadj_sw = x_val; %making the equation zero
                    yadj_sw = y_val;
                else
                    xadj_sw = nudgenode_x(index_n, index_m - 1);
                    yadj_sw = nudgenode_y(index_n, index_m - 1);
                end
                
                %se 
                %index_n + 1, index_m - 1
                
                if index_m == 1 || index_n == 75
                    xadj_se = x_val; %making the equation zero
                    yadj_se = y_val;
                else
                    xadj_se = nudgenode_x(index_n + 1, index_m - 1);
                    yadj_se = nudgenode_y(index_n + 1, index_m - 1);
                end
            end

            
            
            %x
            %%Force Equations enter
            dnw = sqrt((xadj_nw - x_val)^2 + (yadj_nw - y_val)^2);
            dne = sqrt((xadj_ne - x_val)^2 + (yadj_ne - y_val)^2);
            dw = sqrt((xadj_w - x_val)^2 + (yadj_w - y_val)^2);
            de = sqrt((xadj_e - x_val)^2 + (yadj_e - y_val)^2);
            dsw = sqrt((xadj_sw - x_val)^2 + (yadj_sw - y_val)^2);
            dse = sqrt((xadj_se - x_val)^2 + (yadj_se - y_val)^2); 

            %directional k-prime values
            kp_ne = k*(dne - d_o) / dne;
            kp_nw = k*(dnw - d_o) / dnw;
            kp_se = k*(dse - d_o) / dse;
            kp_sw = k*(dsw - d_o) / dsw;
            kp_e = k*(de - d_o) / de;
            kp_w = k*(dw - d_o) / dw;

            
            %Hypothesis 2
            % 
            
            
%Force equation to determine the direction of the nudge
            net_force_x = (kp_ne * (xadj_ne - x_val)) + (kp_nw* (xadj_nw - x_val)) + (kp_se * (xadj_se - x_val)) + (kp_sw * (xadj_sw - x_val)) +  (kp_e * (xadj_e - x_val)) + (kp_w * (xadj_w - x_val));
            net_force_y = (kp_ne * (yadj_ne - y_val)) + (kp_nw* (yadj_nw - y_val)) + (kp_se * (yadj_se - y_val)) + (kp_sw * (yadj_sw - y_val)) +  (kp_e * (yadj_e - y_val)) + (kp_w * (yadj_w - y_val));

          
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
            
%Should there be any unrealistic negative separations of the artery wall,
%eliminate them using the Max function to, where Max(a, b) = a if a > b 
%and Max(a, b) = b otherwise.  Thus for complete collapse the artery wall 
%separation is defined as zero.
%The edge nodes will stay constant
            
            
        end
    end
   
    nudgenode_x = nudge_temp_x;
    nudgenode_y = nudge_temp_y;
    %change index number and add plotting 
    index_force = index_force + 1;
end

%Plots the nudged nodes

x_col = 1;
y_col = 1;
for plotting_2 = 1:18
    hold on
    x_section = nudgenode_x(:,x_col);
    y_section = nudgenode_y(:,y_col);
    figure(f3)
    plot(x_section,y_section,'b.-')
    x_col = x_col + 1;
    y_col = y_col + 1;
end



%This codes for nudge graph
xlabel("Length(cm)");
ylabel("Height(cm)");
title({['Fat test Arm Compression of Tissue with ' ,num2str(cuffwidth),'cm'], ['Blood Pressure Cuff with ',num2str((1 - compressionfactor) * 100),' % Depression']});

hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization for plotting internal artery diameter vs length
