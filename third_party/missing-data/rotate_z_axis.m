function B = rotate_z_axis(A,a)
B = [cos(a), sin(a), 0, 0; -sin(a),cos(a),0,0;0,0,1,0;0,0,0,1]*A;