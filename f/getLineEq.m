function [Line, LineSlopeForm] = getLineEq(P1,P2)
% Find the equation of a Line between 2 points in a plane
% P1 = (x1,y1), P2 = (x2, y2);

nrow1 = size(P1,1);
nrow2 = size(P2,1);
if nrow1 == 1
    P1 = P1';
end
if nrow2 == 1
    P2 = P2';
end


x = P1(1:2,1)';
y = P2(1:2,1)';

Q1 = cart2hom(x);
Q2 = cart2hom(y);

lin = cross(Q1,Q2);
Line = lin/lin(3);

LineSlopeForm = lin/lin(2);

end
