function [vx, vy] = compose_vec_fields(v_second_x,v_second_y,v_first_x,v_first_y)
    
    [x,y] = ndgrid(0:size(v_second_x,1)-1, 0:size(v_second_x,2)-1);
    x_tilda = x + v_first_x;
    y_tilda = y + v_first_y;
    
    comp_x = interpn(x,y,v_second_x,x_tilda,y_tilda,'cubic',0);
    comp_y = interpn(x,y,v_second_y,x_tilda,y_tilda,'cubic',0);
    
    vx = v_first_x + comp_x;
    vy = v_first_y + comp_y;
    
end