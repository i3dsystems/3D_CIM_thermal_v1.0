function [MinT, MaxT] = DrawTran(index, Xmesh, Ymesh, drawT_die, draw, ...
                                 Ta, map, filename, t)
    MinT = min(min(data - 273));
    MaxT = max(max(data - 273));
    
    figure(index);
    range = draw.range;
    if draw.absolutely == 1
        contourf(Xmesh*100, Ymesh*100,drawT_die-273, draw.granularity,'Linestyle','none');
    else
        contourf(Xmesh*100, Ymesh*100,drawT_die-Ta, draw.granularity,'Linestyle','none');
    end
    hold on;
    
    len = size(map, 1);
    
    for i=1:len
        xl = map(i,1);
        width = map(i,3);
        yb = map(i,2);
        height = map(i,4);                    
        rectangle('Position',[xl yb width height]*100);
        hold on;
    end
        
    caxis([range(1), range(2)]);
    axis off;
    axis equal;
    h=colorbar;
    if draw.absolutely == 1
        set(get(h,'Title'),'string','Tjunc','FontSize',16);
    else
        set(get(h,'Title'),'string','Tjunc-Tamb','FontSize',16);
    end
    set(gca,'FontSize',16);
    xlabel('x(cm)');
    ylabel('y(cm)');set(gca,'FontSize',16);

    temp = ['Time: ', num2str(t), 's.'];
    title(temp, 'FontSize', 16);
    frame = getframe(index);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);                
    if t == 0;
        imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.2);
    end
end