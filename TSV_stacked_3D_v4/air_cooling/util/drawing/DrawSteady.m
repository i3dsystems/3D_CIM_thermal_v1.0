function [MaxT, MinT] = DrawSteady(index, Xmesh, Ymesh, drawT_die, draw, ...
                                    Ta, map, name)
    MinT = min(min(drawT_die - 273));
    MaxT = max(max(drawT_die - 273));
    
    figure(index);
    range = draw.range;
    if draw.absolutely == 1
        contourf(Xmesh*100, Ymesh*100,drawT_die-273, draw.granularity,'Linestyle','none');
    else
        contourf(Xmesh*100, Ymesh*100,drawT_die-Ta, draw.granularity,'Linestyle','none');
    end
    hold on;
    
    len = size(map, 1);
%     for i=1:len      
%         if map(i, 6) == 0
%             xl = map(i,1);
%             width = map(i,3);
%             yb = map(i,2);
%             height = map(i,4);                    
%             rectangle('Position',[xl yb width height]*100, 'LineWidth', 1);
%             if isempty(char(name(i))) == 0              
%                 text((xl+width/2)*100, (yb+height/2)*100, char(name(i)), 'HorizontalAlignment','center', 'FontSize', 24)        
%             end
%             hold on;
%         end
%     end
    if draw.clamp == 1
        caxis([range(1), range(2)]);
    end
    axis off;
    axis equal;
    h=colorbar;
    if draw.absolutely == 1
        set(get(h,'Title'),'string','T_{junc}(^{o}C)', 'FontSize',16);
    else
        set(get(h,'Title'),'string','T_{junc}-T_{amb}','FontSize',16);
    end
    set(gca,'FontSize',16);
    xlabel('x(cm)');
    ylabel('y(cm)');set(gca,'FontSize',16);
end