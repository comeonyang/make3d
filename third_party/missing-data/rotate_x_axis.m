function B = rotate_x_axis(A,a)
rot = [1,0,0,0; 0, cos(a), sin(a), 0; 0, -sin(a),cos(a),0;0,0,0,1];
B = rot*A;