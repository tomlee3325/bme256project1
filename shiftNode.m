function [outputArray] = shiftNode(inputArray)
%placeholder will be input for this function
%Iterate x and y positions from placeholder


% Net_Fx > 0 -> shitft = +0.0005
).  If the net x-directed force is > 0 , then nudge the node 5 microns to the right, etc.  If the net x-directed force is < 0 , then nudge the node 5 microns to the left, and similarly for the y-dimension.  
% for possible = 1:500. 
    %per 1 iteration, do we have to test each node's net force?(75*36)
    % per each iteration, we find the positions of all nodes and then
    % calculate the netfore for each node and decide whether they should
    % shift right or left? 
    
%force equation
%Initialization
x_ne = 0;
x_nw = 0 ;
x_e = 0;
x_w = 0;
x_sw = 0;
x_se = 0;



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
updated = []; %this will be an updated array; Reduce naming confusion. 
%iterates through all the nodes
for index_n = 1:75 %x
    for index_m = 1: 18
        %find main node; the node of interest; we will find nw, ne, sw, se
        %based on this node
        %Find node of interest
        x_val = inputArray(index_n, 2*index_m-1);
        y_val = inputArray(index_n, 2*index_m);
        
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
            x_nw = 100;
            x_nw = inputArray(index_n - 1,2*(index_m+1) -1 );
            y_nw = inputArray(index_n - 1, 2* (index_m+1));
            
            %ne
            
        else    %odd
            
            
        %%Force Equations enter
        %net_force_x
        %net_force_y 
                
            
    