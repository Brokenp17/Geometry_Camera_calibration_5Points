function f = getFocalLength(R,T,Points3D, Points2D)

% f = Focal Length
Tx = T(1);
Ty = T(2);
Tz = T(3);
r00 = R(1,1);
r01 = R(1,2);
r02 = R(1,3);
r10 = R(2,1);
r11 = R(2,2);
r12 = R(2,3);
r20 = R(3,1);
r21 = R(3,2);
r22 = R(3,3);

% 3D real world Points = (X Y Z 1)
% 2D image plane Points = (x y 1)

A = [];
B = [];
jj = 1;

for j = 1:length(Points3D)
%     j = j+1
    X = Points3D(1,j);
    Y = Points3D(2,j);
    Z = Points3D(3,j);
    
    x = Points2D(1,j);
    y = Points2D(2,j);
    
    num1 = r00*(X-Tx) + r01*(Y-Ty) + r02*(Z-Tz);
    num2 = r10*(X-Tx) + r11*(Y-Ty) + r12*(Z-Tz);
    den = r20*(X-Tx) + r21*(Y-Ty) + r22*(Z-Tz);
    
    A1 = num1/den;
    A2 = num2/den;
    
    B(jj,:) = x;
    B(jj+1,:) = y;
    A(jj,:) = A1;
    A(jj+1,:) = A2;
    
       jj = jj+2;
end
% LS solution
f = (A.'*A)\(A.'*B);


% end

