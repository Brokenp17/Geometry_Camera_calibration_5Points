function LinePlot(L, color, p00)
%     if isempty(color)
%         color = 'b';
%     end
    
% end
%     if nargin <= 2
%         exsistanceCheck = evalin( 'base', 'exist(''p0'',''var'') == 1' );
%         if exsistanceCheck == 1
%             p00 = evalin('base','p0');
%         end
%     end
%     
    
    
    a = L(1);
    b = L(2);
    c = L(3);
    xlim = get(gca,'XLim');
    ylim = get(gca,'YLim'); 

    if b ~= 0
        m = -a/b;
        k = -c/b;
        FunX = m*xlim + k;
        plot(xlim, FunX, 'LineWidth', 1, 'Color', color);
    else
        x = p00(1);
        line([x x], ylim, 'Color', color);
    end
end

