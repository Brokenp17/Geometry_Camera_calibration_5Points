close all;
clear all;
clc;

addpath f
addpath geom3d %Some function from geom3d Toolbox 
% https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d

Image = imread('File.png');
refObj = imread('NBA_Court.png');
ph = size(refObj,1);

I = Image;
pa = [326.1488; 370.1379];
pb = [731.9750; 388.2113];
pc = [648.1865; 477.2136];
pd = [187.2433; 456.6661];
pe = [225.3878; 165.5465];

pin = [pa pb pc pd]; % 2xN matrix of inputs (Image Points)

% NBA Court Measure in ft
A = [0 50]';
B = [19 50]';
C = [19 33]';
D = [0 33]';
E = [4 28]';
Ze = 13; % Top Backboard from ground in ft
pout = [A B C D]; % 2xN matrix of output (Court Model Points)

% Fit homography from image point to ground plane
H1 = fitgeotrans(pout(1:2,:)', pin(1:2,:)', 'projective');
H = H1.T';

%% Plot the points on the Image 
figure(1);
imshow(Image);
% axis equal

hold on;
InputPoints = [pa pb pc pd pe];
PlotPoints(InputPoints, 'r', 'o')
text(pa(1),pa(2),' pa', 'Color', 'r');
text(pb(1),pb(2),' pb', 'Color', 'r');
text(pc(1),pc(2),' pc', 'Color', 'r');
text(pd(1),pd(2),' pd', 'Color', 'r');
text(pe(1),pe(2),' pe', 'Color', 'r');

p0 = [size(I,2)/2; size(I,1)/2]; % Principal point coordinate, in middle by Hp
InputPoints = [p0];
PlotPoints(InputPoints, 'c', '+')
text(p0(1),p0(2),' p0', 'Color', 'c');

%% Find the Vanishing Line
hold on;
[h00, h01, h02, h10, h11, h12, h20, h21, h22] = DecomposeMtx(H);

v1 = [h00/h20; h10/h20; 1];
v2 = [h01/h21; h11/h21; 1];

%% Get Theta Angle
Theta = atan((h10*h21-h11*h20)/(h00*h21 - h01*h20))

%%
p00 = [1; p0(2)];
    
% line1 = LineFrom2Points(pa,pd);
% line2 = LineFrom2Points(pb,pc);
% LinePlot(line1, 'y'); 
% LinePlot(line2, 'y');
% v1 = cross(line1,line2);

% line1 = LineFrom2Points(pa,pb);
% line2 = LineFrom2Points(pd,pc);
% LinePlot(line1, 'y');
% LinePlot(line2, 'y');
% v2 = cross(line1,line2);

% v1 = (v1/v1(3));
% v2 = (v2/v2(3));

lvan = cross(v1,v2); % Vanishing Line, Horizon
lvan = lvan/lvan(3);

plot(v1(1),v1(2), 'o', 'MarkerFaceColor', 'g'); text(v1(1),v1(2),' v1', 'Color', 'g')
plot(v2(1),v2(2), 'o', 'MarkerFaceColor', 'g'); text(v2(1),v2(2),' v2', 'Color', 'g')

InputPoints = [v1 v2];
plot([v1(1) v2(1)], [v1(2) v2(2)])

bottomline = LineFrom2Points([1, size(I,1)],[size(I,2), size(I,1)]); % Line at image bottom
mbott = -(bottomline(1)/bottomline(2));
kbott = -(bottomline(3)/bottomline(2));

LinePlot(lvan, 'r', p0);
%  l2 = perpendicular to vanishing line passing through p0 in the image
m1 = -(lvan(1)/lvan(2)); 
k1 = -(lvan(3)/lvan(2));

% l2 = perp line eq through p0 : y = m2*x + k2
m2 = -1/m1; % Slope of l2, perpendicular to vanishing line
k2 = p0(2) - m2*p0(1); % y-intercept of l2

% Equation of the vertical line through p0
lvert = cross([p0; 1],[p0(1); 10; 1]);
lvert = lvert/lvert(3);

pperp = LineIntersect([m1, k1], [m2, k2]);
plot(pperp(1),pperp(2), 'o', 'MarkerFaceColor', 'b'); text(pperp(1),pperp(2),' X', 'Color', 'b');

[l2, l22] = LineFrom2Points(pperp,p0);

% Plot all the Lines and Points remaining
LinePlot(lvan, 'r');
LinePlot(l2, 'r');
LinePlot(lvert, 'b',  p0);

% p1 = intersection between L2 and the bottom of the image
p1 = cross(l2,bottomline);
p1 = p1/p1(3);
plot(p1(1),p1(2), 'o', 'MarkerFaceColor', 'c'); text(p1(1),p1(2),' p1', 'Color', 'c');

%% Project p0 and p1 onto Ground (GROUND PLANE = XY PLANE REAL WORLD)
%{
H is a plane 2 plane mapping, we reproject the image points 
p0 and p1 onto the Ground Plane
%}

P0 = hom2cart((H\[p0;1])')'
% P0 = cart2hom(P0')';

P1 = hom2cart(((H)\p1)')';
% P1 = cart2hom(P1')';

Pe = hom2cart((H\[pe;1])');
% Pe = cart2hom(Pe)';

Eprime = cart2hom([E(1) E(2)])';

%% Compute Rotation Angle Omega
Omega = atan((P1(1)-P0(1))/(P0(2)-P1(2)))

%% Lines on Ground
% [L2, L2slope] = getLineEq(P0,P1);
% m2 = -L2slope(1);
m2 = (P1(2) - P0(2))/(P1(1)-P0(1)); %slope formula

%% Find Position of Camera Centre
% Use 5-th point (the non coplanar one)

E = [E(1); E(2); Ze];

Eprime = [E(1); E(2);0;]; % 3d coordinates

% The line Le connecting Pe and E will pass through camera centre G
% Le = cross(Pe,E);
% Le = Le/Le(3);

% The line L3 formed by connecting Pe and E' will pass 
% through G' = [xG, yG, 0]

% [L3, L3slope] = getLineEq(Pe,Eprime);
% m3 = -L3slope(1);
m3 = (Eprime(2)-Pe(2))/(Eprime(1)-Pe(1));

% Gprime = cross(L2, L3);
% Gprime = Gprime/Gprime(3);
% Gprime = Gprime';

% Tx_ = Gprime(1);
Tx = (m2*P0(1) - m3*Pe(1)- P0(2) + Pe(2))/(m2-m3);

% Ty_ = Gprime(2);
Ty = m2*(Tx - P0(1))+P0(2);
% Ty = m3*(Tx - Pe(1)) + Pe(2)

Gprime = [Tx;Ty;0];

PeGprime = [Pe(1) Pe(2); Gprime(1) Gprime(2)];
d1 = pdist(PeGprime,'euclidean');
PeEprime = [Pe(1) Pe(2); Eprime(1) Eprime(2)];
d2 = pdist(PeEprime,'euclidean');

Tz = Ze*(d1/d2);

G = [Tx; Ty; Tz]; % 3d coordinates

%% Phi angle

Dist = [Tx, Ty; P0(1) P0(2)];
Num = pdist(Dist,'Euclidean')
% Num = sqrt((P0(1)-Tx)^2 + (P0(2)-Ty)^2);
Phi = atan(Num/Tz);
Phi = pi/2-abs(Phi)

%% Rotation Matrix and Translation Vector

R0 = [1   0  0;
      0   0  1;
      0   -1  0];
  
R1 = [cos(Theta)    sin(Theta)  0;
     -sin(Theta)    cos(Theta)  0;
        0               0       1];
R2 = [1     0           0;
      0    cos(Phi)  sin(Phi);
      0   -sin(Phi)  cos(Phi)];

R3 = [cos(Omega)   -sin(Omega)  0;
      sin(Omega)    cos(Omega)  0;
        0               0       1];

R = R1*R2*R3;%*R0;
det(R)

T = [Tx;Ty;Tz]

Points3D = [[A;0;1], [B;0;1], [C;0;1], [D;0;1], [E; 1]];
Points2D = [[pa; 1], [pb; 1], [pc; 1], [pd; 1], [pe; 1]];

f = getFocalLength(R',T,Points3D, Points2D)
f = abs(f);

K = diag([f, f, 1]);
K(1:2,3) = p0;

%% CAMERA MATRIX

P = K*R*[eye(3) -T]

%% Plot Warped Image
% I = flipud(I);
% [Iwarp, ref] = imwarp(I,H1);%,'OutputView',imref2d(size(I)));
% I = flipud(I);
% [Iwarp, ref] = imwarp(I,H1);
% Iwar = flipup(Iwarp);
% figure(); imshow(Iwarp); hold on;

%% Birdsview Model
Pe = Pe';

figure();

hold on;
plot(A(1),A(2), 'o', 'MarkerFaceColor', 'r'); text(A(1),A(2),' A', 'Color', 'r'); hold on;
plot(B(1),B(2), 'o', 'MarkerFaceColor', 'r'); text(B(1),B(2),' B', 'Color', 'r'); hold on;
plot(C(1),C(2), 'o', 'MarkerFaceColor', 'r'); text(C(1),C(2),' C', 'Color', 'r'); hold on;
plot(D(1),D(2), 'o', 'MarkerFaceColor', 'r'); text(D(1),D(2),' D', 'Color', 'r'); hold on;
plot(P0(1),P0(2), 'o', 'MarkerFaceColor', 'b'); text(P0(1),P0(2),' P0', 'Color', 'b'); hold on;
plot(P1(1),P1(2), 'o', 'MarkerFaceColor', 'b'); text(P1(1),P1(2),' P1', 'Color', 'b'); hold on;
plot(Pe(1),Pe(2), 'o', 'MarkerFaceColor', 'c'); text(Pe(1),Pe(2),' Pe', 'Color', 'k'); hold on;
plot(Eprime(1),Eprime(2), 'o', 'MarkerFaceColor', 'r'); text(Eprime(1),Eprime(2),' E"', 'Color', 'r'); hold on;
plot(Gprime(1),Gprime(2), 'o', 'MarkerFaceColor', 'm'); text(Gprime(1),Gprime(2),' G"', 'Color', 'm'); hold on;

InputPoints = [P0(1:2,:) Gprime(1:2,:)];
Plot_Lines_Only(InputPoints, 'b', 'o')
InputPoints = [Pe(1:2,:) Gprime(1:2,:)];
Plot_Lines_Only(InputPoints, 'c', 'o')

% Image Limit Points
I = Image;
Limit1 = hom2cart((H\[1; 1; 1])');
Limit2 = hom2cart((H\[size(I,2); 1; 1])');
Limit3 = hom2cart((H\[size(I,2); size(I,1); 1])');
Limit4 = hom2cart((H\[1; size(I,1); 1])');

InputPoints = [Limit1' Limit4'];
Plot_Lines_and_Points(InputPoints, 'g', 'x')
InputPoints = [Limit2' Limit3'];
Plot_Lines_and_Points(InputPoints, 'g', 'x')
InputPoints = [Limit3' Limit4'];
Plot_Lines_and_Points(InputPoints, 'g', 'x')
InputPoints = [Limit1' Limit2'];
Plot_Lines_and_Points(InputPoints, 'g', 'x')

InputPoints = [Limit3' Gprime(1:2,:)];
Plot_Lines_Only(InputPoints, 'y', 'o')

InputPoints = [Limit4' Gprime(1:2,:)];
Plot_Lines_Only(InputPoints, 'y', 'o')

text(Limit1(1),Limit1(2),' 1', 'Color', 'k'); hold on;
text(Limit2(1),Limit2(2),' 2', 'Color', 'k'); hold on;
text(Limit3(1),Limit3(2),' 3', 'Color', 'k'); hold on;
text(Limit4(1),Limit4(2),' 4', 'Color', 'k'); hold on;

%% 3D PLOT %%
close all;
Pe(3) = 0;
P0(3) = 0;
P1(3) = 0;
Eprime(3) = 0;

figure('color','w');
set(gcf, 'Position', get(0, 'Screensize'));
grid on;
xlabel('X'), ylabel('Y'), zlabel('Z')
hold on;

w = 80;
h = 100;
points=[[0 0 0]' [w 0 0]' [w h 0]' [0 h 0]'];
fill3(points(1,:),points(2,:),points(3,:),'r'); hold on;
alpha(0.1)

A = [A;0]; B=[B;0]; C=[C;0]; D=[D;0];
drawPoint3d(E'); hold on;
text(E(1),E(2), E(3),' E', 'Color', 'k');  hold on;
drawPoint3d(A','marker', '+', 'markerSize', 5); hold on;
text(A(1),A(2),0,' A', 'Color', 'k');  hold on;
drawPoint3d(B'); hold on;
text(B(1),B(2),0,' B', 'Color', 'k');  hold on;
drawPoint3d(C'); hold on;
text(C(1),C(2),0,' C', 'Color', 'k');  hold on;
drawPoint3d(D'); hold on;
text(D(1),D(2),0,' D', 'Color', 'k');  hold on;
drawPoint3d([P0(1),P0(2),P0(3)]); hold on;
text(P0(1),P0(2),P0(3),' P0', 'Color', 'k');  hold on;
drawPoint3d([P1(1),P1(2),P1(3)]); hold on;
text(P1(1),P1(2),P1(3),' P1', 'Color', 'k');  hold on;
drawPoint3d([Gprime(1),Gprime(2) 0]); hold on;
text(Gprime(1),Gprime(2),0,' G"', 'Color', 'k');  hold on;
drawPoint3d([Tx, Ty, Tz]); hold on;
text(Tx, Ty, Tz,' Cam', 'Color', 'r');  hold on;

L2_3d = fitLine3d([P0'; P1']);
L3_3d = fitLine3d([Pe'; Eprime']);
Le_3d = fitLine3d([E'; Pe']);
LE_ed = fitLine3d([E'; Eprime']);
L0_cam = fitLine3d([Tx, Ty, Tz; P0(1), P0(2), 0]);
drawLine3d(L2_3d);  hold on;
drawLine3d(L3_3d);  hold on;
drawLine3d(Le_3d);  hold on;
drawLine3d(LE_ed);  hold on;
drawLine3d(L0_cam);  hold on;

R = R1*R2*R3*R0;
% cam = plotCamera('Location',[Tx Ty Tz],'Orientation',eye(3) ,'Opacity',0, 'size', 5); hold on;
cam = plotCamera('Location',[Tx Ty Tz],'Orientation',R' ,'Opacity',0.1, 'size', 5);
