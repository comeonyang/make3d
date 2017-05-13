function INC = motion_incidence(fo,nframes,npoints)
% Each column is a point.  Rows 2i and 2i+1 are the x and y
%  coordinates of the point in frame i.

INC = ones(2*nframes,npoints);
for i = 1:npoints
  % There's got to be an easy way to do this without looping.
  rn = rand;
  if rn < 2/3
     ol = occlusion_length(fo*nframes,nframes);
     if rn < 1/3
        INC((1+2*nframes-2*ol):(2*nframes),i) = zeros(2*ol,1);
     else
        INC(1:2*ol,i) = zeros(2*ol,1);
     end
  else
    ol1 = occlusion_length(fo*nframes/2,nframes/2);
    ol2 = occlusion_length(fo*nframes/2,nframes/2);
    ol = ol1+ol2;
    INC([1:2*ol1,(1+2*nframes-2*ol2):(2*nframes)],i) = zeros(2*ol,1);
  end
end
