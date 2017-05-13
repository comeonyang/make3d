function R = in_plane_rotation(a)
% generate a 4x4 rotation matrix for in-plane rotation.
R = [cos(a), sin(a), 0, 0; -sin(a),cos(a),0,0;0,0,1,0;0,0,0,1];