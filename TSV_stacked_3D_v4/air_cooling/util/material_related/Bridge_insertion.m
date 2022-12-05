function [Cl, Cv, Cp, De] = Bridge_insertion(Cdt_L, Cdt_V, Cap, Den, system, connect_plane, Material, draw, Layer)
% change the conductivity of meshes which contain TSVs
    fprintf('Bridge insertion: there are %i bridge chips \n', system.bridge.N);
    right_id = system.Nchip;
    layer_id = Layer.Nt; %change it if more package layers are inserted
%     bridge_material = system.layer.bridge(2);
    
        if Layer.cal_type(layer_id) == 3
            k1 = Material.K(system.layer.bridge(2));
            k2 = Material.K(Layer.other(layer_id) );
            cp1 = Material.C(system.layer.bridge(2));
            cp2 = Material.C(Layer.other(layer_id) );
            den1 = Material.D(system.layer.bridge(2)); 
            den2 = Material.D(Layer.other(layer_id) );
            
            k = [k1, k2];
            portion = Layer.portion(layer_id);
            p = [1-portion, portion];
            [cv, cl] = Cal_Con(k, p);
            
            cap = cp1*portion + cp2*(1-portion);
            den = den1*portion + den2*(1-portion);
        else
            cv = Material.K(Base_material(layer_id) );
            cl = Material.K(Base_material(layer_id) );
            cap = Material.C(Base_material(layer_id) );
            den = Material.D(Base_material(layer_id) );
        end
        
%     K = Material.K( bridge_material );
%     C = Material.C( bridge_material );
%     D = Material.D( bridge_material );
    const = connect_plane.Nx*connect_plane.Ny;
    
    if draw == 1
        drawC = ones(connect_plane.Ny, connect_plane.Nx, 2);
        drawC(:,:,1) = reshape(Cdt_L(right_id+1:right_id+const, 1), connect_plane.Nx, connect_plane.Ny)';
        drawC(:,:,2) = reshape(Cdt_V(right_id+1:right_id+const), connect_plane.Nx, connect_plane.Ny)';
    end
    
    for i = 1 : system.bridge.N
        xl = find(abs(connect_plane.Xmesh-system.bridge.xl(i))<1e-12);
        xr = find(abs(connect_plane.Xmesh-system.bridge.xr(i))<1e-12);
        yb = find(abs(connect_plane.Ymesh-system.bridge.yb(i))<1e-12);
        yt = find(abs(connect_plane.Ymesh-system.bridge.yt(i))<1e-12);

        Cdt = reshape(Cdt_L(right_id+1:right_id+const, 1), connect_plane.Nx, connect_plane.Ny)';
        Cdt(yb:yt-1, xl:xr-1) = cl;
        Cdt_L(right_id+1:right_id+const, 1) = reshape(Cdt', connect_plane.Nx*connect_plane.Ny, 1);
        Cdt = reshape(Cdt_L(right_id+1:right_id+const, 2), connect_plane.Nx, connect_plane.Ny)';
        Cdt(yb:yt-1, xl:xr-1) = cl;
        Cdt_L(right_id+1:right_id+const, 2) = reshape(Cdt', connect_plane.Nx*connect_plane.Ny, 1);
        if draw == 1
            drawC(yb:yt-1, xl:xr-1, 1) = cl;
        end
        Cdt = reshape(Cdt_V(right_id+1:right_id+const), connect_plane.Nx, connect_plane.Ny)';
        Cdt(yb:yt-1, xl:xr-1) = cv;
        Cdt_V(right_id+1:right_id+const) = reshape(Cdt', connect_plane.Nx*connect_plane.Ny, 1);
        if draw == 1
            drawC(yb:yt-1, xl:xr-1, 2) = cv;
        end
        Cdt = reshape(Cap(right_id+1:right_id+const), connect_plane.Nx, connect_plane.Ny)';
        Cdt(yb:yt-1, xl:xr-1) = cap;
        Cap(right_id+1:right_id+const) = reshape(Cdt', connect_plane.Nx*connect_plane.Ny, 1);

        Cdt = reshape(Den(right_id+1:right_id+const), connect_plane.Nx, connect_plane.Ny)';
        Cdt(yb:yt-1, xl:xr-1) = den;
        Den(right_id+1:right_id+const) = reshape(Cdt', connect_plane.Nx*connect_plane.Ny, 1);    
    end
    
    if draw == 1
        figure(5)
        pcolor(connect_plane.Xmesh*100, connect_plane.Ymesh*100, drawC(:,:,1));
        h=colorbar;
        set(get(h,'Title'),'string','k*m/W','FontSize',16)
        set(gca,'FontSize',16);
        xlabel('x(cm)');
        ylabel('y(cm)');set(gca,'FontSize',16);
        title('plane beneath chip X/Y conductivity');

        figure(6)
        pcolor(connect_plane.Xmesh*100, connect_plane.Ymesh*100, drawC(:,:,2));
        h=colorbar;
        set(get(h,'Title'),'string','k*m/W','FontSize',16)
        set(gca,'FontSize',16);
        xlabel('x(cm)');
        ylabel('y(cm)');set(gca,'FontSize',16);
        title('plane beneath chip Z conductivity');    
    end
    
    fprintf('\n');
    Cl = Cdt_L;
    Cv = Cdt_V;
    Cp = Cap;
    De = Den;
end

