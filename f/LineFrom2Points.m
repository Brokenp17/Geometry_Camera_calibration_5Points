function [Line1, Line2] = LineFrom2Points(A,B)
% Get 2 points A = (x1, y1) and B = (x2, y2) as input and returns the line
% connecting them in aX+bY+c = 0 form and also in slope intercept form as 
% Line2 = [slope, intercept];

nrowA = size(A,1);
nrowB = size(B,1);

if nrowA >= 2
    A = A';
end
if nrowB >= 2
    B = B';
end

x1 = A(1); 
y1 = A(2);
x2 = B(1); 
y2 = B(2);

p1 = [x1; y1; 1];
p2 = [x2; y2; 1];
Line1 = cross(p1,p2);
Line1 = Line1/Line1(3);
%% Check that the line is not parallel to x-axis or y-axis
if (y2 - y1) ~= 0 && (x2 - x1) ~= 0 
    m = (y2 - y1)/(x2 - x1);
    k = y1 - m*x1;
end

%% Line parallel to x-axis
if (y2 - y1) == 0
    m = 0;
    k = y1;   
end
%% Line parallel to y-axis
if (x2 - x1) == 0
    m = x1;
    k = 0;
end
Line2 = [m, k];
end

