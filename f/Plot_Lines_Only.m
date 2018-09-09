function Plot_Lines_Only(A, color, marker)
% input A is a 2xN Matrices;
% Plots the lines from A1 to A2, A2 to A3 etc.
% also plots the line from Aend to A1;

for i = 1:length(A)-1
    j = i+1;
    x1 = A(1,i);
    y1 = A(2,i);
    x2 = A(1,j);
    y2 = A(2,j);

    X = [x1 x2];
    Y = [y1 y2];
    plot(X, Y, color); hold on;
end

% Plot line from last to first element of A
x1 = A(1,1);
y1 = A(2,1);
x2 = A(1,length(A));
y2 = A(2,length(A));

X = [x1 x2];
Y = [y1 y2];
plot(X, Y, color); hold on;

