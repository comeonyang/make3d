function D = randpersp(numpts,numframes)
% produce a random set of projective images.  We'll have no translation, for now.
pts = rand(3,numpts);
trans = rand(3*numframes,3);
pts_trans = trans*pts
D = [];
for i = 1:3:3*numframes
  D = [D; pts_trans(i,:)./pts_trans(i+2,:)];
  D = [D; pts_trans(i+1,:)./pts_trans(i+2,:)];
end
