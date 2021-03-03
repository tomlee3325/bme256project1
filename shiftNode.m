function [updated] = shiftNode(inputArray)
%placeholder will be input for this function
%Iterate x and y positions from placeholder

% do = 0.2 

% Net_Fx > 0 -> shitft = +0.0005
%).  If the net x-directed force is > 0 , then nudge the node 5 microns to the right, etc.  If the net x-directed force is < 0 , then nudge the node 5 microns to the left, and similarly for the y-dimension.  
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

update = [];
update = inputArray;

%not zero its o
d_o = 0.2; %Initial node distance between two nodes after nudge
d = 0.23; %Final node distance

net_force_x = (d .- d_o).*((x_ne - x_value)./d)
%%Net force in x-direction. Zero is equilibrium

index_n = 0; %x-position on the graph
index_m = 0; %y-position on the graph
x_val = 0; %x-position value; temporary storage 
y_val = 0; %y-position value; temporary storage

%indexing x-position(x,2y-1)
%indexing y-position(x,2y)
updated = inputArray; %this will be an updated array; Reduce naming confusion. 
%iterates through all the nodes

for index_n = 1:75 %x
    for index_m = 18: 1 %y
        %find main node; the node of interest; we will find nw, ne, sw, se
        %based on this node
        %Find node of interest
        x_val = inputArray(index_n, 2*index_m-1);
        y_val = inputArray(index_n, 2*index_m);
        if abs(15 - x_val) <= 13.3 
            update(index_n, 2*index_m-1) = y_val * 0.9;
            update(index_n, 2*index_m) = x_val;
        else
            update(index_n, 2*index_m-1) = y_val;
            update(index_n, 2*index_m) = x_val;
        end
    end
end

inputArray()
index_force = 0
while index_force < 500
    
    for index_n = 1:75 %x
        for index_m = 18: 1 %y
            %find main node; the node of interest; we will find nw, ne, sw, se
            %based on this node
            %Find node of interest
            x_val = update(index_n, 2*index_m-1);
            y_val = update(index_n, 2*index_m);

            %use index_m to determine if it is odd row or even row
            %for instance, m represent 1,2,3,4,5,6,7,8,9th row and so on.
            %depending on these row's odd or evenness, the positions of nw, ne
            %and so on

            %find adjacent node positions
            %we will always code in this order
            %nw -> ne -> w-> e -> sw -> se
            if rem(index_m, 2) == 0 %even 
                %nw
                %index_n - 1 , index_m + 1
                x_nw = update((index_n - 1),2*(index_m+1) -1 );
                y_nw = update((index_n - 1), 2* (index_m+1));

                %ne
                %index_n, index_m +1
                x_ne = update(index_n, 2* (index_m + 1) -1);
                x_nw = update(index_n, 2* (index_m + 1));

                %w
                %index_n - 1, index_m
                x_w = update((index_n - 1), 2* index_m - 1);
                y_w = update((index_n-1), 2*index_m); 

                %e
                %index_n + 1, index_m
                x_e = update((index_n + 1), 2* (index_m));
                y_e = update((index_n + 1), 2* (index_m - 1));

                %sw 
                %index_n + 1, index_m - 1
                x_sw = update((index_n + 1), 2* (index_m - 1) - 1);
                y_sw = update((input_n + 1),  2* (index_m - 1)); 

                %se 
                %index_n, index_m - 1
                x_se = update(index_n, 2* (index_m - 1) - 1);
                y_se = update(index_n, 2* (index_m - 1));

            else    %odd
                %nw
                %index_n, index_m + 1
                x_nw = update((index_n),2*(index_m+1) -1 );
                y_nw = update((index_n), 2* (index_m + 1));

                %ne
                %index_n + 1, index_m +1
                x_ne = update(index_n + 1, 2* (index_m + 1) -1);
                y_nw = update(index_n + 1, 2* (index_m + 1));

                %w
                %index_n - 1, index_m
                x_w = update((index_n - 1), 2* index_m - 1);
                y_w = update((index_n-1), 2*index_m); 

                %e
                %index_n + 1, index_m
                x_e = update((index_n + 1), 2* (index_m));
                y_e = update((index_n + 1), 2* (index_m - 1));

                %sw 
                %index_n, index_m - 1
                x_sw = update((index_n), 2* (index_m - 1) - 1);
                y_sw = update((input_n),  2* (index_m - 1)); 

                %se 
                %index_n + 1, index_m - 1
                x_se = update(index_n + 1, 2* (index_m - 1) - 1);
                y_se = update(index_n + 1, 2* (index_m - 1));
            end

            %x
            %%Force Equations enter
            dnw = sqrt((x_nw - x_val)^2 + (y_nw - y_val)^2);
            dne = sqrt((x_ne - x_val)^2 + (y_ne - y_val)^2);
            dw = sqrt((x_w - x_val)^2 + (y_w - y_val)^2);
            de = sqrt((x_e - x_val)^2 + (y_e - y_val)^2);
            dsw = sqrt((x_sw - x_val)^2 + (y_sw - y_val)^2);
            dse = sqrt((x_se - x_val)^2 + (y_se - y_val)^2); 

            d_o = sqrt((x_nw - x_val)^2 + (y_nw - y_val)^2);

%NEED TO UPDATE FORCE EQUATION WITH THE ABOVE VARIABLES
            net_force_x = k((dne - d_o).*((x_ne - x_val)./dne))- k(dnw - d_o).*((x_nw -
            x_val) ./ dnw) + k(dse - d_o)*((x_se-x_val) / dse) - k(dsw-d_o).* ((x_sw -
            x)./ dsw) + k(de - d_o) *((x_e - x_val)./de) - k(dn - do).* ((x_n -
            x_val) ./ dn)

            net_force_y = k((d - d_o).*((y_ne - y_val)./d))- k(d - d_o).*((y_nw -
            y_val) ./ d) + k(d - d_o)*((y_se-y_val) / d) - k(d-d_o).* ((y_sw -
            y)./ d) + k(d - d_o) *((y_e - y_val)./d) - k(d - do).* ((y_n -
            y_val) ./ d)

            if net_force_x > 0 
                updated(index_n, ((index_m * 2)-1)) = x_val + 0.0005;

            elseif net_force_x < 0 
                updated(index_n, ((index_m * 2)-1)) = x_val - 0.0005;

            else
                updated(index_n, ((index_m*2) -1)) = x_val;

            end

            %y
            if net_force_y > 0 
                updated(index_n, (index_m * 2)) = y_val + 0.0005;
            elseif net_force_y < 0 
                updated(index_n, (index_m * 2)) = y_val - 0.0005;
            else
                updated(index_n, (index_m*2)) = y_val;
            end
        end
    end
    
    
    temp_not = update;
    
    %change index number and add plotting 
    index_force = index_force + 1;
end

%k


%         %y
% 
%         end
%     end
%     end
% 
% %update matrix
% %force added
% 
% 
% 
% 
% if net_force_x > 0 
%     updated(index_n, ((index_m * 2)-1)) = x_val + 0.0005;
% elseif net_force_x < 0 
%     updated(index_n, ((index_m * 2)-1)) = x_val - 0.0005;
% else
%     updated(index_n, ((index_m*2) -1)) = x_val;
%     
% end
% 
% 
% if net_force_y > 0 
%     updated(index_n, (index_m * 2)) = y_val + 0.0005;
% elseif net_force_y < 0 
%     updated(index_n, (index_m * 2)) = y_val - 0.0005;
% else
%     updated(index_n, (index_m*2)) = y_val;
%     
% end

    
    
        

    