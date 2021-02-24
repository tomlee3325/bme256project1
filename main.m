
function main(~)

placeholder = zeros(75,36);
delta_x = .2;
delta_y = sqrt(.03);
size_placeholder = size(placeholder);
x_value = 0;
y_value = 0;
nudge_index = 0;
counter_index = 0;
y_index = 0;

for col_index = 1:size_placeholder(2)
    
    if y_index == 2
        y_value = y_value + delta_y;
        y_index = 0;
    end
  
    if counter_index == 4
        x_value = 0;
        counter_index = 0;
        nudge_index = 1;
    elseif nudge_index == 2
        x_value = delta_x / 2;
        %nudge_index = 0;
    end
    
    for row_index = 1:size_placeholder(1)
        if rem(col_index,2) == 0 % even column so y values
            placeholder(row_index,col_index) = y_value;
        else
            placeholder(row_index,col_index) = x_value;
            x_value = delta_x + x_value;
        end
    end
    
    nudge_index = nudge_index + 1;
    counter_index = counter_index + 1;
    y_index = y_index + 1;
end

x_col = 1;
y_col = 2;
for plotting = 1:18
    hold on
    x_section = placeholder(:,x_col);
    y_section = placeholder(:,y_col);
    plot(x_section,y_section,'b.')
    x_col = x_col + 2;
    y_col = y_col + 2;
end
