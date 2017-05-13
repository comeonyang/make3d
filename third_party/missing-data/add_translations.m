function Mt = add_translations(M,ftrans)
% Translation is assumed constant, in a random direction.  ftrans is the total 
% amount.
nframes = size(M,1)/2;
npts = size(M,2);

translation_mag = ftrans/(nframes-1);
% By the last frame, total translation should be ftrans.
translation_dir = 2*pi*rand(1);
translation_vec = translation_mag*[cos(translation_dir), sin(translation_dir)];

A = repeat((0:nframes-1)',npts);
Mt = M;
Mt(1:2:2*nframes,:) = Mt(1:2:2*nframes,:)+translation_vec(1)*A;
Mt(2:2:2*nframes,:) = Mt(2:2:2*nframes,:)+translation_vec(2)*A;
