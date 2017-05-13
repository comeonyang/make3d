function test_colorimetry
% TEST_COLORIMETRY  Test color space conversion functons

[a,b]= meshgrid(0:.01:1,0:.01:1) ;
c = ones(size(a)) ;
d = zeros(size(a)) ;

% 2-D charts
figure(1) ; set(gcf,'Renderer','OpenGL') ; clf ; test('CIE') ;
figure(2) ; set(gcf,'Renderer','OpenGL') ; clf ; test('Adobe') ;

% 3-D charts
figure(3) ; clf ; test3('CIE') ;

% --------------------------------------------------------------------
function test(ws) 
% --------------------------------------------------------------------
% xy + Luv + Lab diagrams
[a,b]= meshgrid(0:.01:1,0:.01:1) ;
c = ones(size(a)) ;
d = zeros(size(a)) ;

draw(a,b,c,ws) ;
draw(c,a,b,ws) ;
draw(b,c,a,ws) ;

% --------------------------------------------------------------------
function test3(ws)
% --------------------------------------------------------------------
% 3-D Lab diagram

[a,b]= meshgrid(0:.01:1,0:.01:1) ;
c = ones(size(a)) ;
d = zeros(size(a)) ;

draw3(a,b,c,ws) ;
draw3(c,a,b,ws) ;
draw3(b,c,a,ws) ;
draw3(a,b,d,ws) ;
draw3(d,a,b,ws) ;
draw3(b,d,a,ws) ;
set(gca,'Box','on') ;

function draw(r,g,b,ws) ;

I = cat(3,r,g,b) ;
J = rgb2xyz(I,ws)  ;

X = J(:,:,1) ;
Y = J(:,:,2) ;
Z = J(:,:,3) ;

x = X./(X+Y+Z) ;
y = Y./(X+Y+Z) ;

J = xyz2luv( cat(3,X,Y,Z)  ) ;
u = J(:,:,2) ;
v = J(:,:,3) ;

J = xyz2lab( cat(3,X,Y,Z)  ) ;
L = J(:,:,1) ;
a = J(:,:,2) ;
b = J(:,:,3) ;

xc = mean(x(:)) ;
yc = mean(y(:)) ;

uc = mean(u(:)) ;
vc = mean(v(:)) ;

ac = mean(a(:)) ;
bc = mean(b(:)) ;

tightsubplot(4,1,'Box','outer') ; hold on ;title([ws ' RGB gamut in xyY']) ;
z0 = zeros(size(X)) ;
h=mesh(x,y,z0) ;
set(h, 'CData', I, 'FaceColor', 'texturemap') ;
campos([xc,yc,1]) ;
camtarget([xc,yc,0]) ;
axis square ;
xlabel('y') ;
ylabel('b') ;

tightsubplot(4,2,'Box','outer') ; hold on ; title([ws ' RGB gamut in Luv']) ;
h=mesh(u,v,z0) ;
set(h, 'CData', I, 'FaceColor', 'texturemap') ;
campos([uc,vc,1]) ;
camtarget([uc,vc,0]) ;
axis square ;
xlabel('u') ;
ylabel('v') ;

tightsubplot(4,3,'Box','outer') ; hold on ; title([ws ' RGB gamut in Lab']) ;
h=mesh(a,b,z0) ;
set(h, 'CData', I, 'FaceColor', 'texturemap') ;
campos([ac,bc,1]) ;
camtarget([ac,bc,0]) ;
axis square ;
xlabel('a') ;
ylabel('b') ;

% --------------------------------------------------------------------
function draw3(r,g,b,ws)
% --------------------------------------------------------------------

I = cat(3,r,g,b) ;
J = rgb2xyz(I,ws)  ;

X = J(:,:,1) ;
Y = J(:,:,2) ;
Z = J(:,:,3) ;

J = xyz2lab( cat(3,X,Y,Z)  ) ;
L = J(:,:,1) ;
a = J(:,:,2) ;
b = J(:,:,3) ;

hold on ; 
title([ws ' RGB gamut in Lab space']) ;
h=mesh(a,b,L) ;
set(h, 'CData', I, 'FaceColor', 'texturemap') ;
xlabel('a') ;
ylabel('b') ;
zlabel('L') ;
axis square ;
grid on ;
