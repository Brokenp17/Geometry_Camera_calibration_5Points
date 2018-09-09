function Plot_Lines_and_Points(A, color, marker)
% input A is a 2xN Matrices;
% Plots the lines from A1 to A2, A2 to A3 etc.
% also plots the line from Aend to A1;

% Plot also the points

for i = 1:length(A)-1
    j = i+1;
    x1 = A(1,i);
    y1 = A(2,i);
    x2 = A(1,j);
    y2 = A(2,j);

    X = [x1 x2];
    Y = [y1 y2];
    plot(X, Y, color); hold on;
%     plot(Limit1(1),Limit1(2), 'o', 'MarkerFaceColor', 'g'); hold on;
    plot(x1, y1, 'Color', color, 'Marker', marker, 'MarkerFaceColor', color); hold on;
end

% Plot line from last to first element of A
x1 = A(1,1);
y1 = A(2,1);
x2 = A(1,length(A));
y2 = A(2,length(A));

X = [x1 x2];
Y = [y1 y2];
plot(X, Y, color); hold on;
plot(x2, y2, 'Color', color, 'Marker', marker, 'MarkerFaceColor', color); hold on;



% plot(Limit1(1),Limit1(2), 'o', 'MarkerFaceColor', 'g'); hold on;
% plot(Limit2(1),Limit2(2), 'o', 'MarkerFaceColor', 'g'); hold on;
% plot(Limit3(1),Limit3(2), 'o', 'MarkerFaceColor', 'g'); hold on;
% plot(Limit4(1),Limit4(2), 'o', 'MarkerFaceColor', 'g'); hold on;


