% TEST_WARP Test image warping routines

% load test image
I = imread('pu.png') ;
I = im2double(I) ;
[M,N,K] = size(I) ;

% dense grid
[x,y] = meshgrid(1:N,1:M);

% control point grid
xr=linspace(1,N,5) ;
yr=linspace(1,M,5) ;
[xc,yc]=meshgrid(xr,yr) ;

% --------------------------------------------------------------------
%                                              test affine defromation
% --------------------------------------------------------------------
A = [1 0 ; 0 1] ;
T = [0;0] ;

figure(100);clf; 
for i=1:10
  % backward projection of dense grid
  [wx,wy]=waffine(A,T,x,y) ;
  J = imwbackward(1:N,1:M,I,wx,wy);
  
  % forward projection of control grid
  [xg,yg]=waffine(inv(A),-inv(A)*T,xc,yc) ;

  % change randomly the affine tf
  A = A + .1*randn(2) ;
  T = T + .1*randn(2,1)  ;
  
  % plot deformed image and grid
  tightsubplot(10,i) ;imagesc(J) ;colormap gray;
  plotgrid(xg,yg,'color','g','linewidth',4) ;axis off;
end

drawnow ;

% --------------------------------------------------------------------
%                                               test thin-plate spline
% --------------------------------------------------------------------
Y= [xc(:)';yc(:)'] ;
U = tps(x,y,Y) ;
Y0=Y ;

figure(101);clf ;
for i=1:10
  % backward projection of dense grid
  [wx,wy] = wtps(U,Y) ;
  J = imwbackward(1:N,1:M,I,wx,wy);

  % forward projection of control grid
  [xg,yg]=witps(xc,yc,Y0,Y) ;

  % change TPS randomly
  Y=Y+3*randn(size(Y)) ;

  % plot deformend image and grid
  tightsubplot(10,i) ;imagesc(J) ;colormap gray; hold on ;
  plotgrid(xg,yg,'color','g','linewidth',4) ; axis off ;
end
