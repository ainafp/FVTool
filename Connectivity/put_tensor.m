function tensor = put_tensor(tensor_base, p_center, p, mode)
    angle = atan((p_center(2)-p(2))/(p_center(1)-p(1)));
    if mode==1
        angle = angle + pi/2;
    elseif mode==2
        angle = 3*pi/2 - angle;
    end
    R = rotate2D(angle);
    tensor = R' * tensor_base * R;
end
