function P = LineIntersect(Line1,Line2)
% Gets as input 2 lines in form : line = [slope, y-intercept] and return
% the crossing point P

m1 = Line1(1);
k1 = Line1(2);
m2 = Line2(1);
k2 = Line2(2);

x = (k2 - k1) / (m1 - m2);
y = m1*x + k1;

P = [x y];

end