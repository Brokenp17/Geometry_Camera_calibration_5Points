function PlotPoints(A, color, marker)
% input A is a 2xN Matrices;
% Plots all the points

for i = 1:size(A,2)
    x1 = A(1,i);
    y1 = A(2,i);
    plot(x1, y1, 'Color', color, 'Marker', marker, 'MarkerFaceColor', color); hold on;
end
