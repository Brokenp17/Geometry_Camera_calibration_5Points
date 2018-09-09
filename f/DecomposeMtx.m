function [a11, a12, a13, a21, a22, a23, a31, a32, a33] = DecomposeMtx(A)
a11 = A(1,1);
a12 = A(1,2);
a13 = A(1,3);

a21 = A(2,1);
a22 = A(2,2);
a23 = A(2,3);

a31 = A(3,1);
a32 = A(3,2);
a33 = A(3,3);

end
