function [Y, C, H] = Matrix_build_FOWLP(chip, system, hc, Cdt_L, Cdt_V, Cap, Den, ...
                               Layer,row, column, value, var)
    t_start = tic;
    fprintf('build Stiffness Matrix\n');
    Thick_layer = Layer.thick;

    h_s = hc.up;
    h_mf = hc.down;
    h_d = hc.d;
    h = hc.side;
    Ta = hc.Ta;
       
    pointer = 1;      

    tic;
    
    layer_id = 1;
        
    pack_xmesh = system.hsp.Xmesh;
    pack_ymesh = system.hsp.Ymesh;

    gridNx_pack = system.hsp.Nx;
    gridNy_pack = system.hsp.Ny;
    
    %% heat spreader top layer
    checker = 1;
    t_check=tic;
    for i = 1:1:gridNx_pack
        for j = 1:1:gridNy_pack
            lz_d = 0;
            lz_u = Thick_layer(layer_id);
            lz = (lz_u+lz_d)/2;

            if(i>1) 
                x1 = pack_xmesh(i)-pack_xmesh(i-1);
            else
                x1 = 0;
            end

            if(i<gridNx_pack)
                x2 = pack_xmesh(i+1) - pack_xmesh(i);
            else
                x2 = 0;
            end

            % y direction
            if(j>1) 
                y1 = pack_ymesh(j)-pack_ymesh(j-1);
            else
                y1 = 0;
            end

            if(j<gridNy_pack)
                y2 = pack_ymesh(j+1) - pack_ymesh(j);
            else
                y2 = 0;
            end

            lx = (x1+x2)/2;
            ly = (y1+y2)/2;

            const = 0;
            id = i+(j-1)*gridNx_pack+const;
            w = id-1;
            e = id+1;
            s = id - gridNx_pack;
            n = id + gridNx_pack;

            t = id + gridNx_pack*gridNy_pack;
            t_s = t - gridNx_pack;                         
            
            if( ( i==1 || i==gridNx_pack ) && ( j==1 || j==gridNy_pack ) )
                if i == 1 && j == 1
                        ke_u = Cdt_L(id, 1);kn_u = Cdt_L(id, 2);kz_u = Cdt_V(id);                    

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u*y2*lz_u/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u*x2*lz_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x2*y2/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = Den(id)*Cap(id)*lz_u*lx*ly/2;
                        value(pointer, 3) = 1/4*(h*y2*lz_u+h*x2*lz_u+h_s*x2*y2);
                        pointer=pointer+1;                    

                elseif i == gridNx_pack && j == 1
                        kw_u = Cdt_L(id-1, 1);kn_u = Cdt_L(id-1, 2);kz_u = Cdt_V(id-1);

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u*y2*lz_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u*x1*lz_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x1*y2/lz_u;         
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = Den(id-1)*Cap(id-1)*lz_u*lx*ly/2;
                        value(pointer, 3) = 1/4*(h*y2*lz_u+h*x1*lz_u+h_s*x1*y2);
                        pointer=pointer+1;                      

                elseif i == 1 && j == gridNy_pack
                        ke_u = Cdt_L(s, 1);ks_u = Cdt_L(s, 2);kz_u = Cdt_V(s);

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u*y1*lz_u/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u*x2*lz_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x2*y1/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = Den(s)*Cap(s)*lz_u*lx*ly/2;
                        value(pointer, 3) = 1/4*(h*y1*lz_u+h*x2*lz_u+h_s*x2*y1);
                        pointer=pointer+1;                     

                else
                        kw_u = Cdt_L(s-1, 1);ks_u = Cdt_L(s-1, 2);kz_u = Cdt_V(s-1);

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*(kw_u*y1*lz_u/x1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*(ks_u*x1*lz_u/y1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*(kz_u*x1*y1/lz_u);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = Den(s-1)*Cap(s-1)*lz_u*lx*ly/2;
                        value(pointer, 3) = 1/4*(h*y1*lz_u+h*x1*lz_u+h_s*x1*y1);
                        pointer=pointer+1;                        

                end
            elseif i ~= 1 && i~= gridNx_pack && j ~= 1 && j~=gridNy_pack
                kw_u = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_u/2;
                ks_u = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_u/2;
                ke_u = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_u/2;
                kn_u = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_u/2;            
                kz_u = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;
                
                row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(2*x2);
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = 0 - value(pointer, 1);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(2*y2);
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = temp - value(pointer, 1);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(2*x1);
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = temp - value(pointer, 1);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(2*y1);
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = temp - value(pointer, 1);
                pointer=pointer+1;           

                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = temp - value(pointer, 1);
                pointer=pointer+1;           

                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                value(pointer, 2) = 1/8 *lz_u * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                            Den(s-1)*Cap(s-1)*x1*y1 + ...
                                            Den(s)*Cap(s)*x2*y1 + ...
                                            Den(id)*Cap(id)*x2*y2);
                value(pointer, 3) = h_s*lx*ly;
                pointer=pointer+1;  

            else
                if i == 1
                        ke_u = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_u/2;
                        kn_u = Cdt_L(id, 2)*x2*lz_u;ks_u = Cdt_L(s, 2)*x2*lz_u;
                        kz_u = (Cdt_V(id)*x2*y2+Cdt_V(s)*x2*y1)/2;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/(x2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(s)*Cap(s)*x2*y1 + ...
                                                    Den(id)*Cap(id)*x2*y2);
                        value(pointer, 3) = h*ly*lz+h_s*lx*ly;
                        pointer=pointer+1;                      

                elseif i == gridNx_pack
                        kw_u = (Cdt_L(id-1, 1)*y2+Cdt_L(s-1, 1)*y1)*lz_u/2;
                        kn_u = Cdt_L(id-1, 2)*x1*lz_u;ks_u = Cdt_L(s-1, 2)*x1*lz_u;
                        kz_u = (Cdt_V(id-1)*x1*y2+Cdt_V(s-1)*x1*y1)/2;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                    Den(s-1)*Cap(s-1)*x1*y1);
                        value(pointer, 3) = h*ly*lz+h_s*lx*ly;
                        pointer=pointer+1;                       

                elseif j == gridNy_pack
                        ks_u = (Cdt_L(s, 2)*x2+Cdt_L(s-1, 2)*x1)*lz_u/2;
                        kw_u = Cdt_L(s-1, 1)*y1*lz_u;ke_u = Cdt_L(s, 1)*y1*lz_u;
                        kz_u = (Cdt_V(s)*x2*y1+Cdt_V(s-1)*x1*y1)/2;

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                    Den(s)*Cap(s)*x2*y1);
                        value(pointer, 3) = h*lx*lz + h_s*lx*ly;
                        pointer=pointer+1;

                else
                        kn_u = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_u/2;
                        kw_u = Cdt_L(id-1, 1)*y2*lz_u;ke_u = Cdt_L(id, 1)*y2*lz_u;
                        kz_u = (Cdt_V(id)*x2*y2+Cdt_V(id-1)*x1*y2)/2; 

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                    Den(id)*Cap(id)*x2*y2);
                        value(pointer, 3) = h*lx*lz+h_s*lx*ly;
                        pointer=pointer+1;                      

                end
            end
        end
    end
    fprintf('heat spreader top done: %f\n', toc(t_check));
    %% other heat spreader layers
    t_check=tic;
    for layer_id = 2:1:Layer.hsp
        checker = checker + 1;
        for i = 1:1:gridNx_pack
            for j = 1:1:gridNy_pack
                lz_u = Thick_layer(layer_id-1);       
                lz_d = Thick_layer(layer_id);
                lz = (lz_u+lz_d)/2;

                % x direction
                if(i>1) 
                    x1 = pack_xmesh(i)-pack_xmesh(i-1);
                else
                    x1 = 0;
                end

                if(i<gridNx_pack)
                    x2 = pack_xmesh(i+1) - pack_xmesh(i);
                else
                    x2 = 0;
                end

                % y direction
                if(j>1) 
                    y1 = pack_ymesh(j)-pack_ymesh(j-1);
                else
                    y1 = 0;
                end

                if(j<gridNy_pack)
                    y2 = pack_ymesh(j+1) - pack_ymesh(j);
                else
                    y2 = 0;
                end

                lx = (x1+x2)/2;
                ly = (y1+y2)/2;

                const = gridNx_pack*gridNy_pack*(layer_id-1);
                id = i+(j-1)*gridNx_pack+const;
                w = id-1;
                e = id+1;
                s = id - gridNx_pack;
                n = id + gridNx_pack;                
                b = id + gridNx_pack*gridNy_pack;
                
                t = id - gridNx_pack*gridNy_pack;
                t_s = t - gridNx_pack;

                if( ( i==1 || i==gridNx_pack ) && ( j==1 || j==gridNy_pack ))
                    if i == 1 && j == 1
                            ke_d = Cdt_L(id, 1)*y2*lz_d;kn_d = Cdt_L(id, 2)*x2*lz_d;kz_d = Cdt_V(id)*x2*y2;  
                            ke_u = Cdt_L(t, 1)*y2*lz_u;kn_u = Cdt_L(t, 2)*x2*lz_u;kz_u = Cdt_V(t)*x2*y2;

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * Den(id)*Cap(id)*x2*y2 + ...
                                                1/8 *lz_u * Den(t)*Cap(t)*x2*y2;
                            value(pointer, 3) = 1/2*(h*y2*lz+h*x2*lz);
                            pointer=pointer+1;                    

                    elseif i == gridNx_pack && j == 1
                            kw_d = Cdt_L(id-1, 1)*y2*lz_d;kn_d = Cdt_L(id-1, 2)*x1*lz_d;kz_d = Cdt_V(id-1)*x1*y2;
                            kw_u = Cdt_L(t-1, 1)*y2*lz_u;kn_u = Cdt_L(t-1, 2)*x1*lz_u;kz_u = Cdt_V(t-1)*x1*y2;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                1/8 *lz_u * Den(t-1)*Cap(t-1)*x1*y2;
                            value(pointer, 3) = 1/2*(h*y2*lz + h*x1*lz);
                            pointer=pointer+1;                      

                    elseif i == 1 && j == gridNy_pack
                            ke_d = Cdt_L(s, 1)*y1*lz_d;ks_d = Cdt_L(s, 2)*x2*lz_d;kz_d = Cdt_V(s)*y1*x2;
                            ke_u = Cdt_L(t_s, 1)*y1*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;kz_u = Cdt_V(t_s)*y1*x2;

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * Den(s)*Cap(s)*x2*y1 + ...
                                                1/8 *lz_u * Den(t_s)*Cap(t_s)*x2*y1; 
                            value(pointer, 3) = 1/2*(h*y1*lz+h*x2*lz);
                            pointer=pointer+1;  

                    else
                            kw_d = Cdt_L(s-1, 1)*y1*lz_d;ks_d = Cdt_L(s-1, 2)*x1*lz_d;kz_d = Cdt_V(s-1)*x1*y1; 
                            kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;kz_u = Cdt_V(t_s-1)*x1*y1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;                    

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                1/8 *lz_u * Den(t_s-1)*Cap(t_s-1)*x1*y1;
                            value(pointer, 3) = 1/2*(h*y1*lz+h*x1*lz);
                            pointer=pointer+1;                      
                    end

                elseif i ~= 1 && i~= gridNx_pack && j ~= 1 && j~=gridNy_pack
                    kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                    ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                    ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                    kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                    kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                    kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                    ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                    ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                    kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;                       
                    kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;                                                             

                    row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/x2-1/2*ke_d/x2;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = 0 - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2-1/2*kn_d/y2;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1-1/2*kw_d/x1;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1-1/2*ks_d/y1;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                    value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                        Den(s-1)*Cap(s-1)*x1*y1 + ...
                                        Den(s)*Cap(s)*x2*y1 + ...
                                        Den(id)*Cap(id)*x2*y2);
                    value(pointer, 2) = 1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                        Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                        Den(t_s)*Cap(t_s)*x2*y1 + ...
                                        Den(t)*Cap(t)*x2*y2) + value(pointer, 2);
                    value(pointer, 3) = 0;
                    pointer=pointer+1;                      

                else
                    if i == 1
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = Cdt_L(id, 2)*x2*lz_d;ks_d = Cdt_L(s, 2)*x2*lz_d;
                            kz_d = (Cdt_V(id)*y2+Cdt_V(s)*y1)*x2/2;

                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = Cdt_L(t, 2)*x2*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;
                            kz_u = (Cdt_V(t)*y2+Cdt_V(t_s)*y1)*x2/2;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/x2-1/2*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = h*ly*lz;
                            pointer=pointer+1;                    

                    elseif i == gridNx_pack
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            kn_d = Cdt_L(id-1, 2)*x1*lz_d;ks_d = Cdt_L(s-1, 2)*x1*lz_d;
                            kz_d = (Cdt_V(id-1)*y2+Cdt_V(s-1)*y1)*x1/2;  

                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            kn_u = Cdt_L(t-1, 2)*x1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;
                            kz_u = (Cdt_V(t-1)*y2+Cdt_V(t_s-1)*y1)*x1/2;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1-1/2*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                Den(s-1)*Cap(s-1)*x1*y1) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                Den(t_s-1)*Cap(t_s-1)*x1*y1);
                            value(pointer, 3) = h*ly*lz;
                            pointer=pointer+1;                      

                    elseif j == gridNy_pack
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            kw_d = Cdt_L(s-1, 1)*y1*lz_d;ke_d = Cdt_L(s, 1)*y1*lz_d;
                            kz_d = (Cdt_V(s)*x2+Cdt_V(s-1)*x1)*y1/2;  

                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ke_u = Cdt_L(t_s, 1)*y1*lz_u;
                            kz_u = (Cdt_V(t_s)*x2+Cdt_V(t_s-1)*x1)*y1/2;                         

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1-1/2*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) =  1/8 *lz_d * ( Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                 Den(s)*Cap(s)*x2*y1) + ...
                                                 1/8 *lz_u * ( Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                 Den(t_s)*Cap(t_s)*x2*y1);
                            value(pointer, 3) = h*lx*lz;
                            pointer=pointer+1; 

                    else
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;
                            kw_d = Cdt_L(id-1, 1)*y2*lz_d;ke_d = Cdt_L(id, 1)*y2*lz_d;
                            kz_d = (Cdt_V(id)*x2+Cdt_V(id-1)*x1)*y2/2;  

                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;
                            kw_u = Cdt_L(t-1, 1)*y2*lz_u;ke_u = Cdt_L(t, 1)*y2*lz_u;
                            kz_u = (Cdt_V(t)*x2+Cdt_V(t-1)*x1)*y2/2;                        

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2-1/2*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = h*lx*lz;
                            pointer=pointer+1;                     
                    end
                end
            end
        end  
    end
    fprintf('other heat spreader layers done: %f\n', toc(t_check));
    layer_id = Layer.hsp+1;

    %use the following 'for loop' when you have multiple chips placed 
    %laterally (such as a 2.5D config).
%     box.hsp = zeros(system.chip.N, 4);
%     box.pack = zeros(system.chip.N, 4);
%     for k = 1 : 1 : system.chip.N               
%         box.hsp(k, :) = ChipInPlane(chip(k), system.hsp);
%         box.pack(k, :) = ChipInPlane(chip(k), system.pack);
%     end
    
    box.hsp = ChipInPlane(chip(system.chip.N), system.hsp);
    box.pack = ChipInPlane(chip(system.chip.N), system.pack);
    
    %% heat spreader bottom layer
    checker = checker + 1;
    t_check=tic;
    for i = 1:1:gridNx_pack
        for j=1:1:gridNy_pack
            lz_u = Thick_layer(layer_id-1);           
            if(i>1) 
                x1 = pack_xmesh(i)-pack_xmesh(i-1);
            else
                x1 = 0;
            end

            if(i<gridNx_pack)
                x2 = pack_xmesh(i+1) - pack_xmesh(i);
            else
                x2 = 0;
            end

            % y direction
            if(j>1) 
                y1 = pack_ymesh(j)-pack_ymesh(j-1);
            else
                y1 = 0;
            end

            if(j<gridNy_pack)
                y2 = pack_ymesh(j+1) - pack_ymesh(j);
            else
                y2 = 0;
            end

            lx = (x1+x2)/2;
            ly = (y1+y2)/2;
            const = gridNx_pack*gridNy_pack*Layer.hsp;
            id = i+(j-1)*gridNx_pack+const;
            w = id-1;
            e = id+1;
            s = id - gridNx_pack;
            n = id + gridNx_pack;
            t = id - gridNx_pack*gridNy_pack;
            t_s = t - gridNx_pack;
            
            if (id == 282)
                q = 1;
            end
            
%             [b, lz_d] = hsp_id2pack_id(Layer, system.Nhsp, system.chip, i, j, box.hsp);
%             xl = box.hsp(1); xr = box.hsp(2); yb = box.hsp(3); yt = box.hsp(4);
            
            die_id = box2chip (i, j, box.hsp);
            if die_id ~= 0
                [b, lz_d] = hsp_id2chip_id(die_id, Layer, system.Nhsp, chip, i, j, box.hsp(die_id,:));
                temp = box.hsp(die_id,:);
                xl = temp(1); xr = temp(2); yb = temp(3); yt = temp(4);
            end
            
            if( ( i==1 || i==gridNx_pack ) && ( j==1 || j==gridNy_pack ))
                if i == 1 && j == 1
                        ke_u = Cdt_L(t, 1);kn_u = Cdt_L(t, 2);kz_u = Cdt_V(t);
                        
                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u*y2*lz_u/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u*x2*lz_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x2*y2/lz_u;                        
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * Den(t)*Cap(t)*x2*y2; 
                        value(pointer, 3) = 1/4*(h*y2*x2+h*x2*lz_u+h*y2*lz_u);
                        pointer=pointer+1;

                elseif i == gridNx_pack && j == 1
                        kw_u = Cdt_L(t-1, 1);kn_u = Cdt_L(t-1, 2);kz_u = Cdt_V(t-1);

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u*y2*lz_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u*x1*lz_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x1*y2/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;                 

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2);
                        value(pointer, 3) = 1/4*(h*y2*lz_u+h*x1*lz_u+h*x1*y2);
                        pointer=pointer+1;                     

                elseif i == 1 && j == gridNy_pack
                        ke_u = Cdt_L(t_s, 1);ks_u = Cdt_L(t_s, 2);kz_u = Cdt_V(t_s); 

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u*y1*lz_u/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u*x2*lz_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x2*y1/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * Den(t_s)*Cap(t_s)*x2*y1; 
                        value(pointer, 3) = 1/4*(h*y1*lz_u+h*x2*lz_u+h*x2*y1);
                        pointer=pointer+1;                       

                else
                        kw_u = Cdt_L(t_s-1, 1);ks_u = Cdt_L(t_s-1, 2);kz_u = Cdt_V(t_s-1);

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u*y1*lz_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u*x1*lz_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x1*y1/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * (Den(t_s-1)*Cap(t_s-1)*x1*y1); 
                        value(pointer, 3) = 1/4*(h*y1*lz_u+h*x1*lz_u+h*x1*y1);
                        pointer=pointer+1;                                           
                end
            elseif i ~= 1 && i ~= gridNx_pack && j ~= 1 && j ~= gridNy_pack
%                 if(  i >= xl && i <= xr && j >= yb && j <= yt )
                if(  die_id ~= 0 )
                    if (i==xl || i==xr) && (j==yb || j==yt)
                        if i==xl && j==yb
                            ke_d = Cdt_L(id, 1)*y2*lz_d;kn_d = Cdt_L(id, 2)*x2*lz_d;kz_d = Cdt_V(id)*x2*y2;
                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                            kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_d/(4*x2)-ke_u/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_d/(4*y2)-kn_u/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/(lz_u);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/(4*lz_d);            
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)= temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 0;
                            pointer=pointer+1;
                        elseif i==xl && j==yt
                            ke_d = Cdt_L(s, 1)*y1*lz_d;ks_d = Cdt_L(s, 2)*x2*lz_d;kz_d = Cdt_V(s)*x2*y1;                        
                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                            kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;
                                                
                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_d/(4*x2)-ke_u/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_d/(4*y1)-ks_u/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/(4*lz_d);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;             
                            value(pointer, 2) = 1/8 *lz_d * ( Den(s)*Cap(s)*x2*y1) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 0;
                            pointer=pointer+1;

                        elseif i==xr && j==yb
                            kw_d = Cdt_L(id-1, 1)*y2*lz_d;kn_d = Cdt_L(id-1, 2)*x1*lz_d;kz_d = Cdt_V(id-1)*x1*y2;
                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                            kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;                                                  

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_d/(4*y2)-kn_u/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_d/(4*x1)-kw_u/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/(4*lz_d);            
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 0;
                            pointer=pointer+1;
                        else
                            kw_d = Cdt_L(s-1, 1)*y1*lz_d;ks_d = Cdt_L(s-1, 1)*x1*lz_d;kz_d = Cdt_V(s-1)*x1*y1;
                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                            kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;
                                                                          
                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_d/(4*x1)-kw_u/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_d/(4*y1)-ks_u/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/(4*lz_d);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)= temp; 
                            value(pointer, 2) = 1/8 *lz_d * ( Den(s-1)*Cap(s-1)*x1*y1) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 0;
                            pointer=pointer+1;
                        end

                    elseif i~=xl && i~=xr && j~=yb && j~=yt
                        kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                        ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                        ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                        kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                        kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                        kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                        ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                        ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                        kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                        kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;
                                               
                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/x2-1/2*ke_d/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2-1/2*kn_d/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1-1/2*kw_d/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1-1/2*ks_d/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                    Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                    Den(s)*Cap(s)*x2*y1 + ...
                                                    Den(id)*Cap(id)*x2*y2) + ...
                                            1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                    Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                    Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                    Den(t)*Cap(t)*x2*y2);
                        value(pointer, 3) = 0;
                        pointer=pointer+1;                          
                    else
                        if j == yb
                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                            kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;

                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;
                            kw_d = Cdt_L(id-1, 1)*y2*lz_d;ke_d = Cdt_L(id, 1)*y2*lz_d;
                            kz_d = (Cdt_V(id)*x2*y2+Cdt_V(id-1)*x1*y2)/2;

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_d/(4*x2)-ke_u/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_d/(2*y2)-kn_u/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_d/(4*x1)-kw_u/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/(2*lz_d);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 0;
                            pointer=pointer+1;                  

                        elseif j == yt
                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                            kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;

                            ks_d = (Cdt_L(s, 2)*x2+Cdt_L(s-1, 2)*x1)*lz_d/2;
                            kw_d = Cdt_L(s-1, 1)*y1*lz_d;ke_d = Cdt_L(s, 1)*y1*lz_d;
                            kz_d = (Cdt_V(s)*x2*y1+Cdt_V(s-1)*x1*y1)/2;                                                     

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_d/(4*x2)-ke_u/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_d/(4*x1)-kw_u/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1); 
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_d/(2*y1)-ks_u/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/(2*lz_d);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 0;
                            pointer=pointer+1;                  

                        elseif i == xl
                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                            kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;

                            ke_d = (Cdt_L(id, 1)*y2+Cdt_L(s, 1)*y1)*lz_d/2;
                            kn_d = Cdt_L(id, 2)*x2*lz_d;ks_d = Cdt_L(s, 2)*x2*lz_d;
                            kz_d = (Cdt_V(id)*x2*y2+Cdt_V(s)*x2*y1)/2;                                                   
                                                       
                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(2*x2)-ke_d/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_d/(4*y2)-kn_u/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_d/(4*y1)-ks_u/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/(2*lz_d);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 0;
                            pointer=pointer+1;                                           

                        else
                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                            kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;

                            kw_d = (Cdt_L(id-1, 1)*y2+Cdt_L(s-1, 1)*y1)*lz_d/2;
                            kn_d = Cdt_L(id-1, 2)*x1*lz_d;ks_d = Cdt_L(s-1, 2)*x1*lz_d;
                            kz_d = (Cdt_V(id-1)*x1*y2+Cdt_V(s-1)*x1*y1)/2;                                                   

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_d/(4*y2)-kn_u/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_d/(2*x1)-kw_u/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_d/(4*y1)-ks_u/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/(2*lz_d); 
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;    
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 0;
                            pointer=pointer+1;  
                        end
                    end
                else
                    kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                    ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                    ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                    kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                    kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;                    

                    row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(2*x2);
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = 0 - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(2*y2);
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(2*x1);
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(2*y1);
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;           

                    row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;           

                    row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;   
                    value(pointer, 2) = 1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                Den(t)*Cap(t)*x2*y2);
                    value(pointer, 3) = h*lx*ly;
                    pointer=pointer+1;                
                end        
            else
                if i == 1
                        ke_u = (Cdt_L(t, 1)*y2+Cdt_L(t_s, 1)*y1)*lz_u/2;
                        kn_u = Cdt_L(t, 2)*x2*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;
                        kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t_s)*x2*y1)/2;                        

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/(2*y1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/(2*y2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/(x2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;   
                        value(pointer, 2) = 1/8 *lz_u * ( Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                    Den(t)*Cap(t)*x2*y2); 
                        value(pointer, 3) = 1/2*(h*ly*x2+h*ly*lz_u);
                        pointer=pointer+1;                    

                elseif i == gridNx_pack
                        kw_u = (Cdt_L(t-1, 1)*y2+Cdt_L(t_s-1, 1)*y1)*lz_u/2;
                        kn_u = Cdt_L(t-1, 2)*x1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;
                        kz_u = (Cdt_V(t-1)*x1*y2+Cdt_V(t_s-1)*x1*y1)/2;                                                                                         

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/(2*y1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/(x1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/(2*y2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                    Den(t_s-1)*Cap(t_s-1)*x1*y1);                          
                        value(pointer, 3) = 1/2*(h*ly*x1+h*ly*lz_u);
                        pointer=pointer+1;                          

                elseif j == gridNy_pack
                        ks_u = (Cdt_L(t_s, 2)*x2+Cdt_L(t_s-1, 2)*x1)*lz_u/2;
                        kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ke_u = Cdt_L(t_s, 1)*y1*lz_u;
                        kz_u = (Cdt_V(t_s)*x2*y1+Cdt_V(t_s-1)*x1*y1)/2;                

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/(2*x2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/(2*x1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;                   

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;    
                        value(pointer, 2) = 1/8 *lz_u * ( Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                    Den(t_s)*Cap(t_s)*x2*y1);
                        value(pointer, 3) = 1/2*(h*lx*lz_u+h*lx*y1);
                        pointer=pointer+1;                                              

                else
                        kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;
                        kw_u = Cdt_L(t-1, 1)*y2*lz_u;ke_u = Cdt_L(t, 1)*y2*lz_u;
                        kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t-1)*x1*y2)/2;                

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/(2*x2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/(2*x1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                    Den(t)*Cap(t)*x2*y2); 
                        value(pointer, 3) = 1/2*(h*lx*lz_u+h*lx*y2);
                        pointer=pointer+1;                                          
                end
            end       
        end
    end    
    fprintf('heat spreader bot done: %f\n', toc(t_check));
    %% chip layers
    
%     gridNx_pack = system.pack.Nx;
%     gridNy_pack = system.pack.Ny; 
%     pack_xmesh = system.pack.Xmesh;
%     pack_ymesh = system.pack.Ymesh;    
    
    t_check=tic;
%     for k = 1 : 1 : system.chip.N       
        chip_xmesh = chip(1).Xmesh;
        chip_ymesh = chip(1).Ymesh;    
        gridNx_chip = chip(1).Nx;
        gridNy_chip = chip(1).Ny;
%         start = sum(Layer.N(1:k-1))-(k-1)+Layer.hsp+1; %need to change this to a simpler way
%         
%         index_start = 0;
%         for l= 1 : 1 : k-1
%             index_start = index_start + chip(l).Nx*chip(l).Ny*(Layer.N(l)-1);
%         end
%         
%         for layer_id = start+1:1:Layer.N(k)+start-1
        for layer_id = Layer.hsp+2 : 1 : (Layer.Nt - Layer.pack)
            checker = checker + 1;
            for i = 1:1:gridNx_chip
                for j = 1:1:gridNy_chip
%                     lz_u = Thick_layer(layer_id-1+k-1);       
%                     lz_d = Thick_layer(layer_id+k-1);
                    lz_u = Thick_layer(layer_id-1);
                    lz_d = Thick_layer(layer_id);
                    lz = (lz_u+lz_d)/2;

                    if(i>1) 
                        x1 = chip_xmesh(i)-chip_xmesh(i-1);
                    else
                        x1 = 0;
                    end

                    if(i<gridNx_chip)
                        x2 = chip_xmesh(i+1) - chip_xmesh(i);
                    else
                        x2 = 0;
                    end

                    % y direction
                    if(j>1) 
                        y1 = chip_ymesh(j)-chip_ymesh(j-1);
                    else
                        y1 = 0;
                    end

                    if(j<gridNy_chip)
                        y2 = chip_ymesh(j+1) - chip_ymesh(j);
                    else
                        y2 = 0;
                    end

                    lx = (x1+x2)/2;
                    ly = (y1+y2)/2;                

                    const = gridNx_chip*gridNy_chip*(layer_id-(Layer.hsp+2)) + system.Nhsp;
%                     const = gridNx_chip*gridNy_chip*(layer_id-1-start) + system.Nhsp + index_start;
                    id = i+(j-1)*gridNx_chip+const;
                    w = id-1;
                    e = id+1;
                    s = id - gridNx_chip;
                    n = id + gridNx_chip;
                    
%                     b = id+gridNx_chip*gridNy_chip;

                    if layer_id == Layer.hsp+2
                        temp = box.hsp;
                        xl = temp(1); xr = temp(2); yb = temp(3); yt = temp(4);                       
                        t = (i-1+xl)+(j-2+yb)*system.hsp.Nx+system.hsp.Nx*system.hsp.Ny*Layer.hsp;
                        t_s = t - system.hsp.Nx;
                    else                        
                        t = id-gridNx_chip*gridNy_chip;
                        t_s = t - gridNx_chip;
                    end
                    
                    if (t == 367)
                        q = 1;
                    end
                    
                    if layer_id == Layer.Nt-Layer.pack
                        temp = box.pack;
                        xl = temp(1); xr = temp(2); yb = temp(3); yt = temp(4);                        
                        b = (i-1+xl)+(j-2+yb)*system.pack.Nx + system.Nchip;
                    else
                        b = id+gridNx_chip*gridNy_chip;
                    end
                    

%                     if layer_id == Layer.N(k)+start-1   
%                         temp = box.pack;
%                         xl = temp(1); xr = temp(2); yb = temp(3); yt = temp(4);                        
%                         b = (i-1+xl)+(j-2+yb)*system.pack.Nx + system.Nchip;
%                     else
%                         b = id+gridNx_chip*gridNy_chip;
%                     end
%                     
%                     if layer_id == start + 1
%                         temp = box.hsp;
%                         xl = temp(1); xr = temp(2); yb = temp(3); yt = temp(4);                       
%                         t = (i-1+xl)+(j-2+yb)*system.hsp.Nx+system.hsp.Nx*system.hsp.Ny*Layer.hsp;
%                         t_s = t - system.hsp.Nx;
%                     else                        
%                         t = id-gridNx_chip*gridNy_chip;
%                         t_s = t - gridNx_chip;
%                     end
                    
                    if( ( i==1 || i==gridNx_chip ) && ( j==1 || j==gridNy_chip ))
                        if i == 1 && j == 1
                                ke_d = Cdt_L(id, 1)*y2*lz_d;kn_d = Cdt_L(id, 2)*x2*lz_d;kz_d = Cdt_V(id)*x2*y2;  
                                ke_u = Cdt_L(t, 1)*y2*lz_u;kn_u = Cdt_L(t, 2)*x2*lz_u;kz_u = Cdt_V(t)*x2*y2;

                                row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = 0 - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;           

                                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                                value(pointer, 2) =  1/8 *lz_d * Den(id)*Cap(id)*x2*y2 + ...
                                                     1/8 *lz_u * Den(t)*Cap(t)*x2*y2; 
                                value(pointer, 3) = 1/2*(h*y2*lz+h*x2*lz);
                                pointer=pointer+1;                    

                        elseif i == gridNx_chip && j == 1
                                kw_d = Cdt_L(id-1, 1)*y2*lz_d;kn_d = Cdt_L(id-1, 2)*x1*lz_d;kz_d = Cdt_V(id-1)*x1*y2;
                                kw_u = Cdt_L(t-1, 1)*y2*lz_u;kn_u = Cdt_L(t-1, 2)*x1*lz_u;kz_u = Cdt_V(t-1)*x1*y2;

                                row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = 0 - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                                value(pointer, 2) = 1/8 *lz_d * Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                    1/8 *lz_u * Den(t-1)*Cap(t-1)*x1*y2;  
                                value(pointer, 3) = 1/2*(h*y2*lz+h*x1*lz);
                                pointer=pointer+1;                      

                        elseif i == 1 && j == gridNy_chip
                                ke_d = Cdt_L(s, 1)*y1*lz_d;ks_d = Cdt_L(s, 2)*x2*lz_d;kz_d = Cdt_V(s)*y1*x2;
                                ke_u = Cdt_L(t_s, 1)*y1*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;kz_u = Cdt_V(t_s)*y1*x2;

                                row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = 0 - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;            

                                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                                value(pointer, 2) = 1/8 *lz_d * Den(s)*Cap(s)*x2*y1 + ...
                                                    1/8 *lz_u * Den(t_s)*Cap(t_s)*x2*y1;
                                value(pointer, 3) = 1/2*(h*y1*lz+h*x2*lz);
                                pointer=pointer+1;

                        else
                                kw_d = Cdt_L(s-1, 1)*y1*lz_d;ks_d = Cdt_L(s-1, 2)*x1*lz_d;kz_d = Cdt_V(s-1)*x1*y1; 
                                kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;kz_u = Cdt_V(t_s-1)*x1*y1;

                                row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = 0 - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;                    

                                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                                value(pointer, 2) = 1/8 *lz_d * Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                    1/8 *lz_u * Den(t_s-1)*Cap(t_s-1)*x1*y1;
                                value(pointer, 3) = 1/2*(h*y1*lz+h*x1*lz);
                                pointer=pointer+1;                      
                        end

                    elseif i ~= 1 && i~= gridNx_chip && j ~= 1 && j~=gridNy_chip
                        kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                        ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                        ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                        kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                        kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                        kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                        ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                        ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                        kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                        kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/x2-1/2*ke_d/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2-1/2*kn_d/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1-1/2*kw_d/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1-1/2*ks_d/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                    Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                    Den(s)*Cap(s)*x2*y1 + ...
                                                    Den(id)*Cap(id)*x2*y2) + ...
                                            1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                    Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                    Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                    Den(t)*Cap(t)*x2*y2);        
                        value(pointer, 3) = 0;
                        pointer=pointer+1;                      

                    else
                        if i == 1
                                ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                                kn_d = Cdt_L(id, 2)*x2*lz_d;ks_d = Cdt_L(s, 2)*x2*lz_d;
                                kz_d = (Cdt_V(id)*y2+Cdt_V(s)*y1)*x2/2;

                                ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                                kn_u = Cdt_L(t, 2)*x2*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;
                                kz_u = (Cdt_V(t)*y2+Cdt_V(t_s)*y1)*x2/2;                        

                                row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = 0 - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/x2-1/2*ke_d/x2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;            

                                row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                                value(pointer, 2) = 1/8 *lz_d * ( Den(s)*Cap(s)*x2*y1 + ...
                                                            Den(id)*Cap(id)*x2*y2) + ...
                                                    1/8 *lz_u * ( Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                            Den(t)*Cap(t)*x2*y2);        

                                value(pointer, 3) = h*ly*lz;
                                pointer=pointer+1;                    

                        elseif i == gridNx_chip
                                kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                                kn_d = Cdt_L(id-1, 2)*x1*lz_d;ks_d = Cdt_L(s-1, 2)*x1*lz_d;
                                kz_d = (Cdt_V(id-1)*y2+Cdt_V(s-1)*y1)*x1/2;  

                                kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                                kn_u = Cdt_L(t-1, 2)*x1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;
                                kz_u = (Cdt_V(t-1)*y2+Cdt_V(t_s-1)*y1)*x1/2;                         

                                row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = 0 - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1-1/2*kw_d/x1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;            

                                row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;            

                                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                                value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                            Den(s-1)*Cap(s-1)*x1*y1) + ...
                                                    1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                            Den(t_s-1)*Cap(t_s-1)*x1*y1);        
                                value(pointer, 3) = h*ly*lz;
                                pointer=pointer+1;                      

                        elseif j == gridNy_chip
                                ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                                kw_d = Cdt_L(s-1, 1)*y1*lz_d;ke_d = Cdt_L(s, 1)*y1*lz_d;
                                kz_d = (Cdt_V(s)*x2+Cdt_V(s-1)*x1)*y1/2;  

                                ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                                kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ke_u = Cdt_L(t_s, 1)*y1*lz_u;
                                kz_u = (Cdt_V(t_s)*x2+Cdt_V(t_s-1)*x1)*y1/2;                         

                                row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = 0 - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1-1/2*ks_d/y1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;            

                                row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                                value(pointer, 2) = 1/8 *lz_d * ( Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                            Den(s)*Cap(s)*x2*y1) + ...
                                                    1/8 *lz_u * ( Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                            Den(t_s)*Cap(t_s)*x2*y1);        
                                value(pointer, 3) = h*lx*lz;
                                pointer=pointer+1; 

                        else
                                kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;
                                kw_d = Cdt_L(id-1, 1)*y2*lz_d;ke_d = Cdt_L(id, 1)*y2*lz_d;
                                kz_d = (Cdt_V(id)*x2+Cdt_V(id-1)*x1)*y2/2;  

                                kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;
                                kw_u = Cdt_L(t-1, 1)*y2*lz_u;ke_u = Cdt_L(t, 1)*y2*lz_u;
                                kz_u = (Cdt_V(t)*x2+Cdt_V(t-1)*x1)*y2/2;                        

                                row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = 0 - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2-1/2*kn_d/y2;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;            

                                row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;

                                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u; 
                                value(pointer, 2) = 0; value(pointer, 3) = 0;
                                temp = temp - value(pointer, 1);
                                pointer=pointer+1;            

                                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                                value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                            Den(id)*Cap(id)*x2*y2) + ...
                                                    1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                            Den(t)*Cap(t)*x2*y2);        
                                value(pointer, 3) = h*lx*lz;
                                pointer=pointer+1;                     

                        end
                    end
                end
            end
        end
%     end
    fprintf('other chip layers done: %f\n', toc(t_check));
    
    pack_xmesh = system.pack.Xmesh;
    pack_ymesh = system.pack.Ymesh;

    gridNx_pack = system.pack.Nx;
    gridNy_pack = system.pack.Ny;


    %% top of package layer
    %because all the TIM layers are merged to one
    layer_id = Layer.Nt-Layer.pack+1;
    checker = checker + 1;
    t_check=tic;
    for i = 1:1:gridNx_pack
        for j=1:1:gridNy_pack
            lz_d = Thick_layer(layer_id);
            
            if(i>1) 
                x1 = pack_xmesh(i)-pack_xmesh(i-1);
            else
                x1 = 0;
            end

            if(i<gridNx_pack)
                x2 = pack_xmesh(i+1) - pack_xmesh(i);
            else
                x2 = 0;
            end

            % y direction
            if(j>1) 
                y1 = pack_ymesh(j)-pack_ymesh(j-1);
            else
                y1 = 0;
            end

            if(j<gridNy_pack)
                y2 = pack_ymesh(j+1) - pack_ymesh(j);
            else
                y2 = 0;
            end

            lx = (x1+x2)/2;
            ly = (y1+y2)/2;
            const = system.Nchip;
            id = i+(j-1)*gridNx_pack+const;
            w = id-1;
            e = id+1;
            s = id - gridNx_pack;
            n = id + gridNx_pack;
            b = id+gridNx_pack*gridNy_pack;

            die_id = box2chip (i, j, box.pack);
            if die_id ~= 0
                lz_u = Thick_layer(layer_id-1);
                t = pack_id2chip_id(die_id, system.Nchip, chip, i, j, box.pack(die_id,:));
%                 [t, lz_u] = pack_id2chip_id(die_id, Layer, system.Nhsp, chip, i, j, box.pack(die_id,:));
                t_s = t - chip(die_id).Nx;
                temp = box.pack(die_id,:);
                xl = temp(1); xr = temp(2); yb = temp(3); yt = temp(4);
            end            
                        
            if( ( i==1 || i==gridNx_pack ) && ( j==1 || j==gridNy_pack ))
                if i == 1 && j == 1
                        ke_d = Cdt_L(id, 1);kn_d = Cdt_L(id, 2);kz_d = Cdt_V(id);

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_d*y2*lz_d/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_d*x2*lz_d/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d*x2*y2/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_d * (Den(id)*Cap(id)*x2*y2); 
                        value(pointer, 3) = 1/4*(h*y2*lz_d+h*x2*lz_d+h*x2*y2);
                        pointer=pointer+1;

                elseif i == gridNx_pack && j == 1
                        kw_d = Cdt_L(id-1, 1);kn_d = Cdt_L(id-1, 2);kz_d = Cdt_V(id-1);
                        
                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_d*y2*lz_d/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_d*x1*lz_d/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d*x1*y2/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;                 

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2); 
                        value(pointer, 3) = 1/4*(h*y2*lz_d+h*x1*lz_d+h*x1*y2);
                        pointer=pointer+1;                     

                elseif i == 1 && j == gridNy_pack
                        ke_d = Cdt_L(s, 1);ks_d = Cdt_L(s, 2);kz_d = Cdt_V(s);

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_d*y1*lz_d/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_d*x2*lz_d/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d*x2*y1/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_d * ( Den(s)*Cap(s)*x2*y1 );
                        value(pointer, 3) = 1/4*(h*y1*lz_d+h*x2*lz_d+h*x2*y1);
                        pointer=pointer+1;                       

                else
                        kw_d = Cdt_L(s-1, 1);ks_d = Cdt_L(s-1, 2);kz_d = Cdt_V(s-1);

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_d*y1*lz_d/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_d*x1*lz_d/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d*x1*y1/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_d * ( Den(s-1)*Cap(s-1)*x1*y1);
                        value(pointer, 3) = 1/4*(h*y1*lz_d+h*x1*lz_d+h*x1*y1);
                        pointer=pointer+1;                                           

                end
            elseif i ~= 1 && i ~= gridNx_pack && j ~= 1 && j ~= gridNy_pack
                if die_id ~= 0
                    if (i==xl || i==xr) && (j==yb || j==yt)
                        %% four corners
                        if i==xl && j==yb
                            ke_u = Cdt_L(t, 1)*y2*lz_u;kn_u = Cdt_L(t, 2)*x2*lz_u;kz_u = Cdt_V(t)*x2*y2;
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;
                           
                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(4*x2)-ke_d/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(4*y2)-kn_d/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_d/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_d/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/(lz_d);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/(4*lz_u);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)= temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 1/4*(h*y2*lz_u+h*x2*lz_u+h*x1*y2+h*x1*y1+h*x2*y1);
                            pointer=pointer+1;
                            
                        elseif i==xl && j==yt
                            ke_u = Cdt_L(t_s, 1)*y1*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;kz_u = Cdt_V(t_s)*x2*y1;                        
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;                          
                            
                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(4*x2)-ke_d/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_d/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_d/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(4*y1)-ks_d/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/(4*lz_u);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t_s)*Cap(t_s)*x2*y1);
                            value(pointer, 3) = 1/4*(h*y1*lz_u+h*x2*lz_u+h*x1*y2+h*x1*y1+h*x2*y2);
                            pointer=pointer+1;

                        elseif i==xr && j==yb
                            kw_u = Cdt_L(t-1, 1)*y2*lz_u;kn_u = Cdt_L(t-1, 2)*x1*lz_u;kz_u = Cdt_V(t-1)*x1*y2;
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_d/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(4*y2)-kn_d/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(4*x1)-kw_d/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_d/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/(4*lz_u);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp; 
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2);
                            value(pointer, 3) = 1/4*(h*y2*lz_u+h*x1*lz_u+h*x2*y1+h*x1*y1+h*x2*y2);
                            pointer=pointer+1;
                        else
                            kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;kz_u = Cdt_V(t_s-1)*x1*y1;
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;                          
                            
                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_d/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_d/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(4*x1)-kw_d/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(4*y1)-ks_d/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/(4*lz_u);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)= temp;   
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t_s-1)*Cap(t_s-1)*x1*y1);
                            value(pointer, 3) = 1/4*(h*y1*lz_u+h*x1*lz_u+h*x2*y1+h*x1*y2+h*x2*y2);
                            pointer=pointer+1;
                        end
                    %% middle points
                    elseif i~=xl && i~=xr && j~=yb && j~=yt
                        kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                        ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                        ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                        kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                        kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                        kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                        ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                        ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                        kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                        kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/x2-1/2*ke_d/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2-1/2*kn_d/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1-1/2*kw_d/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1-1/2*ks_d/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                    Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                    Den(s)*Cap(s)*x2*y1 + ...
                                                    Den(id)*Cap(id)*x2*y2) + ...
                                            1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                    Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                    Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                    Den(t)*Cap(t)*x2*y2);
                        value(pointer, 3) = 0;
                        pointer=pointer+1;                          
                    else
                        if j == yb
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;
                            kw_u = Cdt_L(t-1, 1)*y2*lz_u;ke_u = Cdt_L(t, 1)*y2*lz_u;
                            kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t-1)*x1*y2)/2;                                                 
                            
                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(4*x2)-ke_d/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(2*y2)-kn_d/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(4*x1)-kw_d/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_d/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/(2*lz_u);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;      
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t)*Cap(t)*x2*y2);
                            value(pointer, 3) = 1/2*(h*lx*y1+h*lx*lz_u);                            
                            pointer=pointer+1;                  

                        elseif j == yt
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            ks_u = (Cdt_L(t_s, 2)*x2+Cdt_L(t_s-1, 2)*x1)*lz_u/2;
                            kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ke_u = Cdt_L(t_s, 1)*y1*lz_u;
                            kz_u = (Cdt_V(t_s)*x2*y1+Cdt_V(t_s-1)*x1*y1)/2;                        

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(4*x2)-ke_d/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_d/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(4*x1)-kw_d/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1); 
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(2*y1)-ks_d/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/(2*lz_u);  
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2); 
                            value(pointer, 2) = 1/8 *lz_u * ( Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1) + value(pointer, 2);
                            value(pointer, 3) = 1/2*(h*lx*y2+h*lx*lz_u);                             
                            pointer=pointer+1;                  

                        elseif i == xl
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            ke_u = (Cdt_L(t, 1)*y2+Cdt_L(t_s, 1)*y1)*lz_u/2;
                            kn_u = Cdt_L(t, 2)*x2*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;
                            kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t_s)*x2*y1)/2;                        

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_d/(2*x2)-ke_u/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(4*y2)-kn_d/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_d/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(4*y1)-ks_d/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/(2*lz_u); 
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;     
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2); 
                            value(pointer, 2) = 1/8 *lz_u * ( Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2) + value(pointer, 2);
                            value(pointer, 3) = 1/2*(h*ly*x1+h*ly*lz_u);
                            pointer=pointer+1;                                           

                        else
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                            kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                            kw_u = (Cdt_L(t-1, 1)*y2+Cdt_L(t_s-1, 1)*y1)*lz_u/2;
                            kn_u = Cdt_L(t-1, 2)*x1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;
                            kz_u = (Cdt_V(t-1)*x1*y2+Cdt_V(t_s-1)*x1*y1)/2;                                                  
                            
                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_d/(2*x2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(4*y2)-kn_d/(2*y2);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(2*x1)-kw_d/(2*x1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(4*y1)-ks_d/(2*y1);
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/(2*lz_u);     
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;  
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2);
                            value(pointer, 2) = 1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1) + value(pointer, 2);
                            value(pointer, 3) = 1/2*(h*ly*x2+h*ly*lz_u);                            
                            pointer=pointer+1;  
                        end
                    end
                else
                    kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                    ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                    ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                    kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                    kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;                    

                    row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_d/(2*x2);
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = 0 - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_d/(2*y2);
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_d/(2*x1);
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_d/(2*y1);
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;           

                    row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;           

                    row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;   
                    value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                Den(s)*Cap(s)*x2*y1 + ...
                                                Den(id)*Cap(id)*x2*y2);
                    value(pointer, 3) = h*lx*ly;  
                    pointer=pointer+1;                
                end
            else
                if i == 1
                        ke_d = (Cdt_L(id, 1)*y2+Cdt_L(s, 1)*y1)*lz_d/2;
                        kn_d = Cdt_L(id, 2)*x2*lz_d;ks_d = Cdt_L(s, 2)*x2*lz_d;
                        kz_d = (Cdt_V(id)*x2*y2+Cdt_V(s)*x2*y1)/2;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_d/(2*y1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_d/(2*y2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_d/(x2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;    
                        value(pointer, 2) = 1/8 *lz_d * ( Den(s)*Cap(s)*x2*y1 + ...
                                                    Den(id)*Cap(id)*x2*y2);     
                        value(pointer, 3) = 1/2*(h*ly*lz_d+h*x2*ly);
                        pointer=pointer+1;                    

                elseif i == gridNx_pack
                        kw_d = (Cdt_L(id-1, 1)*y2+Cdt_L(s-1, 1)*y1)*lz_d/2;
                        kn_d = Cdt_L(id-1, 2)*x1*lz_d;ks_d = Cdt_L(s-1, 2)*x1*lz_d;
                        kz_d = (Cdt_V(id-1)*x1*y2+Cdt_V(s-1)*x1*y1)/2;                

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_d/(2*y1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_d/(x1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_d/(2*y2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;    
                        value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                    Den(s-1)*Cap(s-1)*x1*y1 );  
                        value(pointer, 3) = 1/2*(h*ly*lz_d+h*x1*ly);
                        pointer=pointer+1;                          

                elseif j == gridNy_pack
                        ks_d = (Cdt_L(s, 2)*x2+Cdt_L(s-1, 2)*x1)*lz_d/2;
                        kw_d = Cdt_L(s-1, 1)*y1*lz_d;ke_d = Cdt_L(s, 1)*y1*lz_d;
                        kz_d = (Cdt_V(s)*x2*y1+Cdt_V(s-1)*x1*y1)/2;                

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_d/(2*x2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_d/(2*x1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_d/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;                   

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;    
                        value(pointer, 2) = 1/8 *lz_d * ( Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                    Den(s)*Cap(s)*x2*y1); 
                        value(pointer, 3) = 1/2*(h*lx*lz_d+h*lx*y1);
                        pointer=pointer+1;                                              

                else
                        kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;
                        kw_d = Cdt_L(id-1, 1)*y2*lz_d;ke_d = Cdt_L(id, 1)*y2*lz_d;
                        kz_d = (Cdt_V(id)*x2*y2+Cdt_V(id-1)*x1*y2)/2;                

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_d/(2*x2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_d/(2*x1);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_d/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;     

                        row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                    Den(id)*Cap(id)*x2*y2); 
                        value(pointer, 3) = 1/2*(h*lx*lz_d+h*lx*y2);
                        pointer=pointer+1;                                          
                end
            end       
        end
    end
    fprintf('package top layer done: %f\n', toc(t_check));
    %% the other bottom layers
    t_check=tic;
    for layer_id = Layer.Nt-Layer.pack+2 : Layer.Nt
        checker = checker + 1;
        for i = 1:1:gridNx_pack
            for j = 1:1:gridNy_pack
                lz_u = Thick_layer(layer_id-1);       
                lz_d = Thick_layer(layer_id);
                lz = (lz_u+lz_d)/2;

                % x direction
                if(i>1) 
                    x1 = pack_xmesh(i)-pack_xmesh(i-1);
                else
                    x1 = 0;
                end

                if(i<gridNx_pack)
                    x2 = pack_xmesh(i+1) - pack_xmesh(i);
                else
                    x2 = 0;
                end

                % y direction
                if(j>1) 
                    y1 = pack_ymesh(j)-pack_ymesh(j-1);
                else
                    y1 = 0;
                end

                if(j<gridNy_pack)
                    y2 = pack_ymesh(j+1) - pack_ymesh(j);
                else
                    y2 = 0;
                end

                lx = (x1+x2)/2;
                ly = (y1+y2)/2;

                const = system.Nchip + gridNx_pack*gridNy_pack*(layer_id-(Layer.Nt+1-Layer.pack));
%                 const = system.Nhsp+gridNx_pack*gridNy_pack*(layer_id - (Layer.hsp+1) -1);
                id = i+(j-1)*gridNx_pack+const;
                w = id-1;
                e = id+1;
                s = id - gridNx_pack;
                n = id + gridNx_pack;                
                b = id + gridNx_pack*gridNy_pack;
                
                t = id - gridNx_pack*gridNy_pack;
                t_s = t - gridNx_pack;

                if( ( i==1 || i==gridNx_pack ) && ( j==1 || j==gridNy_pack ))
                    if i == 1 && j == 1
                            ke_d = Cdt_L(id, 1)*y2*lz_d;kn_d = Cdt_L(id, 2)*x2*lz_d;kz_d = Cdt_V(id)*x2*y2;  
                            ke_u = Cdt_L(t, 1)*y2*lz_u;kn_u = Cdt_L(t, 2)*x2*lz_u;kz_u = Cdt_V(t)*x2*y2;

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;           

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) =  1/8 *lz_d * Den(id)*Cap(id)*x2*y2 + ...
                                                 1/8 *lz_u * Den(t)*Cap(t)*x2*y2; 
                            value(pointer, 3) = 1/2*(h*y2*lz+h*x2*lz);
                            pointer=pointer+1;                    

                    elseif i == gridNx_pack && j == 1
                            kw_d = Cdt_L(id-1, 1)*y2*lz_d;kn_d = Cdt_L(id-1, 2)*x1*lz_d;kz_d = Cdt_V(id-1)*x1*y2;
                            kw_u = Cdt_L(t-1, 1)*y2*lz_u;kn_u = Cdt_L(t-1, 2)*x1*lz_u;kz_u = Cdt_V(t-1)*x1*y2;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                1/8 *lz_u * Den(t-1)*Cap(t-1)*x1*y2;  
                            value(pointer, 3) = 1/2*(h*y2*lz+h*x1*lz);
                            pointer=pointer+1;                      

                    elseif i == 1 && j == gridNy_pack
                            ke_d = Cdt_L(s, 1)*y1*lz_d;ks_d = Cdt_L(s, 2)*x2*lz_d;kz_d = Cdt_V(s)*y1*x2;
                            ke_u = Cdt_L(t_s, 1)*y1*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;kz_u = Cdt_V(t_s)*y1*x2;

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * Den(s)*Cap(s)*x2*y1 + ...
                                                1/8 *lz_u * Den(t_s)*Cap(t_s)*x2*y1;
                            value(pointer, 3) = 1/2*(h*y1*lz+h*x2*lz);
                            pointer=pointer+1;

                    else
                            kw_d = Cdt_L(s-1, 1)*y1*lz_d;ks_d = Cdt_L(s-1, 2)*x1*lz_d;kz_d = Cdt_V(s-1)*x1*y1; 
                            kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;kz_u = Cdt_V(t_s-1)*x1*y1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/4*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;                    

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                1/8 *lz_u * Den(t_s-1)*Cap(t_s-1)*x1*y1;
                            value(pointer, 3) = 1/2*(h*y1*lz+h*x1*lz);
                            pointer=pointer+1;                      
                    end

                elseif i ~= 1 && i~= gridNx_pack && j ~= 1 && j~=gridNy_pack
                    kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                    ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                    ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                    kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;            
                    kz_d = (Cdt_V(s-1)*x1*y1+Cdt_V(id-1)*x1*y2+Cdt_V(s)*x2*y1+Cdt_V(id)*x2*y2)/4;

                    kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                    ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                    ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                    kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                    kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;
                   
                    row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/x2-1/2*ke_d/x2;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = 0 - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2-1/2*kn_d/y2;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1-1/2*kw_d/x1;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1-1/2*ks_d/y1;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=b; value(pointer, 1)=-kz_d/lz_d;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                    value(pointer, 2) = 0; value(pointer, 3) = 0;
                    temp = temp - value(pointer, 1);
                    pointer=pointer+1;

                    row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                    value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                Den(s)*Cap(s)*x2*y1 + ...
                                                Den(id)*Cap(id)*x2*y2) + ...
                                        1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                Den(t)*Cap(t)*x2*y2);        
                    value(pointer, 3) = 0;
                    pointer=pointer+1;                      

                else
                    if i == 1
                            ke_d = (Cdt_L(s, 1)*y1+Cdt_L(id, 1)*y2)*lz_d/2;
                            kn_d = Cdt_L(id, 2)*x2*lz_d;ks_d = Cdt_L(s, 2)*x2*lz_d;
                            kz_d = (Cdt_V(id)*y2+Cdt_V(s)*y1)*x2/2;

                            ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                            kn_u = Cdt_L(t, 2)*x2*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;
                            kz_u = (Cdt_V(t)*y2+Cdt_V(t_s)*y1)*x2/2;                        

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/x2-1/2*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(s)*Cap(s)*x2*y1 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                        Den(t)*Cap(t)*x2*y2);        

                            value(pointer, 3) = h*ly*lz;
                            pointer=pointer+1;                    

                    elseif i == gridNx_pack
                            kw_d = (Cdt_L(s-1, 1)*y1+Cdt_L(id-1, 1)*y2)*lz_d/2;
                            kn_d = Cdt_L(id-1, 2)*x1*lz_d;ks_d = Cdt_L(s-1, 2)*x1*lz_d;
                            kz_d = (Cdt_V(id-1)*y2+Cdt_V(s-1)*y1)*x1/2;  

                            kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                            kn_u = Cdt_L(t-1, 2)*x1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;
                            kz_u = (Cdt_V(t-1)*y2+Cdt_V(t_s-1)*y1)*x1/2;                         

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1-1/4*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1-1/2*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2-1/4*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(s-1)*Cap(s-1)*x1*y1) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t_s-1)*Cap(t_s-1)*x1*y1);        
                            value(pointer, 3) = h*ly*lz;
                            pointer=pointer+1;                      

                    elseif j == gridNy_pack
                            ks_d = (Cdt_L(s-1, 2)*x1+Cdt_L(s, 2)*x2)*lz_d/2;
                            kw_d = Cdt_L(s-1, 1)*y1*lz_d;ke_d = Cdt_L(s, 1)*y1*lz_d;
                            kz_d = (Cdt_V(s)*x2+Cdt_V(s-1)*x1)*y1/2;  

                            ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                            kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ke_u = Cdt_L(t_s, 1)*y1*lz_u;
                            kz_u = (Cdt_V(t_s)*x2+Cdt_V(t_s-1)*x1)*y1/2;                         

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1-1/2*ks_d/y1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(s-1)*Cap(s-1)*x1*y1 + ...
                                                        Den(s)*Cap(s)*x2*y1) + ...
                                                1/8 *lz_u * ( Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                        Den(t_s)*Cap(t_s)*x2*y1);        
                            value(pointer, 3) = h*lx*lz;
                            pointer=pointer+1; 

                    else
                            kn_d = (Cdt_L(id, 2)*x2+Cdt_L(id-1, 2)*x1)*lz_d/2;
                            kw_d = Cdt_L(id-1, 1)*y2*lz_d;ke_d = Cdt_L(id, 1)*y2*lz_d;
                            kz_d = (Cdt_V(id)*x2+Cdt_V(id-1)*x1)*y2/2;  

                            kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;
                            kw_u = Cdt_L(t-1, 1)*y2*lz_u;ke_u = Cdt_L(t, 1)*y2*lz_u;
                            kz_u = (Cdt_V(t)*x2+Cdt_V(t-1)*x1)*y2/2;                        

                            row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2-1/4*ke_d/x2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = 0 - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1-1/4*kw_d/x1;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2-1/2*kn_d/y2;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=b; value(pointer, 1)=-1/2*kz_d/lz_d;
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;

                            row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u; 
                            value(pointer, 2) = 0; value(pointer, 3) = 0;
                            temp = temp - value(pointer, 1);
                            pointer=pointer+1;            

                            row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                            value(pointer, 2) = 1/8 *lz_d * ( Den(id-1)*Cap(id-1)*x1*y2 + ...
                                                        Den(id)*Cap(id)*x2*y2) + ...
                                                1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                        Den(t)*Cap(t)*x2*y2);        
                            value(pointer, 3) = h*lx*lz;
                            pointer=pointer+1;                     

                    end
                end
            end
        end  
    end
    fprintf('other package layers done: %f\n', toc(t_check));

    
    
    %% bottom of the chip
    checker = checker + 1;
    layer_id = Layer.Nt+1;
    t_check=tic;
    for i = 1:1:gridNx_pack
        for j = 1:1:gridNy_pack
            lz_d = 0;
            lz_u = Thick_layer(layer_id-1);
            lz = (lz_u+lz_d)/2;

            if(i>1) 
                x1 = pack_xmesh(i)-pack_xmesh(i-1);
            else
                x1 = 0;
            end

            if(i<gridNx_pack)
                x2 = pack_xmesh(i+1) - pack_xmesh(i);
            else
                x2 = 0;
            end

            % y direction
            if(j>1) 
                y1 = pack_ymesh(j)-pack_ymesh(j-1);
            else
                y1 = 0;
            end

            if(j<gridNy_pack)
                y2 = pack_ymesh(j+1) - pack_ymesh(j);
            else
                y2 = 0;
            end

            lx = (x1+x2)/2;
            ly = (y1+y2)/2;

            const = system.Nchip+gridNx_pack*gridNy_pack*(Layer.pack);
            id = i+(j-1)*gridNx_pack+const;        
            w = id-1;
            e = id+1;
            s = id - gridNx_pack;
            n = id + gridNx_pack;

            t = id-gridNx_pack*gridNy_pack;
            t_s = t - gridNx_pack; 

            if( ( i==1 || i==gridNx_pack ) && ( j==1 || j==gridNy_pack ) )
                if i == 1 && j == 1
                        ke_u = Cdt_L(t, 1);kn_u = Cdt_L(t, 2);kz_u = Cdt_V(t);

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u*y2*lz_u/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u*x2*lz_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x2*y2/lz_u;         
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * Den(t)*Cap(t)*x2*y2;
                        value(pointer, 3) = 1/4*(h*y2*lz_u+h*x2*lz_u+h_d*x2*y2);
                        pointer=pointer+1;                    

                elseif i == gridNx_pack && j == 1
                        kw_u = Cdt_L(t-1, 1);kn_u = Cdt_L(t-1, 2);kz_u = Cdt_V(t-1);

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u*y2*lz_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u*x1*lz_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x1*y2/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * Den(t-1)*Cap(t-1)*x1*y2;    
                        value(pointer, 3) = 1/4*(h*y2*lz_u + h*x1*lz_u + h_d*x1*y2);
                        pointer=pointer+1;                      

                elseif i == 1 && j == gridNy_pack
                        ke_u = Cdt_L(t_s, 1);ks_u = Cdt_L(t_s, 2);kz_u = Cdt_V(t_s);

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u*y1*lz_u/x2;
                        value(pointer, 2) = 0; 
                        value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u*x2*lz_u/y1;
                        value(pointer, 2) = 0; 
                        value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*kz_u*x2*y1/lz_u;     
                        value(pointer, 2) = 0; 
                        value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * Den(t_s)*Cap(t_s)*x2*y1; 
                        value(pointer, 3) = 1/4*(h*y1*lz_u+h*x2*lz_u+h_d*x2*y1);
                        pointer=pointer+1;                     

                else
                        kw_u = Cdt_L(t_s-1, 1);ks_u = Cdt_L(t_s-1, 2);kz_u = Cdt_V(t_s-1);

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*(kw_u*y1*lz_u/x1);
                        value(pointer, 2) = 0; 
                        value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*(ks_u*x1*lz_u/y1);
                        value(pointer, 2) = 0; 
                        value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/4*(kz_u*x1*y1/lz_u);
                        value(pointer, 2) = 0; 
                        value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;          
                        value(pointer, 2) = 1/8 *lz_u * Den(t_s-1)*Cap(t_s-1)*x1*y1;  
                        value(pointer, 3) = 1/4*(h*y1*lz_u+h*x1*lz_u+h_d*x1*y1);
                        pointer=pointer+1;                        

                end
            elseif i ~= 1 && i~= gridNx_pack && j ~= 1 && j~=gridNy_pack
                kw_u = (Cdt_L(t_s-1, 1)*y1+Cdt_L(t-1, 1)*y2)*lz_u/2;
                ks_u = (Cdt_L(t_s-1, 2)*x1+Cdt_L(t_s, 2)*x2)*lz_u/2;
                ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;            
                kz_u = (Cdt_V(t_s-1)*x1*y1+Cdt_V(t-1)*x1*y2+Cdt_V(t_s)*x2*y1+Cdt_V(t)*x2*y2)/4;                

                row(pointer)=id; column(pointer)=e; value(pointer, 1)=-ke_u/(2*x2);
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = 0 - value(pointer, 1);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=n; value(pointer, 1)=-kn_u/(2*y2);
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = temp - value(pointer, 1);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=w; value(pointer, 1)=-kw_u/(2*x1);
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = temp - value(pointer, 1);
                pointer=pointer+1;

                row(pointer)=id; column(pointer)=s; value(pointer, 1)=-ks_u/(2*y1);
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = temp - value(pointer, 1);
                pointer=pointer+1;           

                row(pointer)=id; column(pointer)=t; value(pointer, 1)=-kz_u/lz_u;
                value(pointer, 2) = 0; value(pointer, 3) = 0;
                temp = temp - value(pointer, 1);
                pointer=pointer+1;           

                row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                value(pointer, 2) = 1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                            Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                            Den(t_s)*Cap(t_s)*x2*y1 + ...
                                            Den(t)*Cap(t)*x2*y2); 
                if i >= xl && i<= xr && j>=yb && j<=yt  
                    value(pointer, 3) = h_mf*lx*ly;
                else
                    value(pointer, 3) = h_d*lx*ly;
                end
                pointer=pointer+1;  

            else
                if i == 1
                        ke_u = (Cdt_L(t_s, 1)*y1+Cdt_L(t, 1)*y2)*lz_u/2;
                        kn_u = Cdt_L(t, 2)*x2*lz_u;ks_u = Cdt_L(t_s, 2)*x2*lz_u;
                        kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t_s)*x2*y1)/2;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/2*ke_u/(x2);
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(t_s)*Cap(t_s)*x2*y1 + ...
                                                    Den(t)*Cap(t)*x2*y2);  
                        value(pointer, 3) = h*ly*lz+h_d*lx*ly;
                        pointer=pointer+1;                      

                elseif i == gridNx_pack
                        kw_u = (Cdt_L(t-1, 1)*y2+Cdt_L(t_s-1, 1)*y1)*lz_u/2;
                        kn_u = Cdt_L(t-1, 2)*x1*lz_u;ks_u = Cdt_L(t_s-1, 2)*x1*lz_u;
                        kz_u = (Cdt_V(t-1)*x1*y2+Cdt_V(t_s-1)*x1*y1)/2;                     

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/4*ks_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/2*kw_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/4*kn_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                          Den(t_s-1)*Cap(t_s-1)*x1*y1);
                        value(pointer, 3) = h*ly*lz+h_d*lx*ly;
                        pointer=pointer+1;                       

                elseif j == gridNy_pack
                        ks_u = (Cdt_L(t_s, 2)*x2+Cdt_L(t_s-1, 2)*x1)*lz_u/2;
                        kw_u = Cdt_L(t_s-1, 1)*y1*lz_u;ke_u = Cdt_L(t_s, 1)*y1*lz_u;
                        kz_u = (Cdt_V(t_s)*x2*y1+Cdt_V(t_s-1)*x1*y1)/2;

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=s; value(pointer, 1)=-1/2*ks_u/y1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) = 1/8 *lz_u * ( Den(t_s-1)*Cap(t_s-1)*x1*y1 + ...
                                                          Den(t_s)*Cap(t_s)*x2*y1);  
                        value(pointer, 3) = h*lx*lz+h_d*lx*ly;
                        pointer=pointer+1;                        

                else
                        kn_u = (Cdt_L(t, 2)*x2+Cdt_L(t-1, 2)*x1)*lz_u/2;
                        kw_u = Cdt_L(t-1, 1)*y2*lz_u;ke_u = Cdt_L(t, 1)*y2*lz_u;
                        kz_u = (Cdt_V(t)*x2*y2+Cdt_V(t-1)*x1*y2)/2;                    

                        row(pointer)=id; column(pointer)=e; value(pointer, 1)=-1/4*ke_u/x2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = 0 - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=w; value(pointer, 1)=-1/4*kw_u/x1;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;

                        row(pointer)=id; column(pointer)=n; value(pointer, 1)=-1/2*kn_u/y2;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;         

                        row(pointer)=id; column(pointer)=t; value(pointer, 1)=-1/2*kz_u/lz_u;
                        value(pointer, 2) = 0; value(pointer, 3) = 0;
                        temp = temp - value(pointer, 1);
                        pointer=pointer+1;           

                        row(pointer)=id; column(pointer)=id; value(pointer, 1)=temp;
                        value(pointer, 2) =  1/8 *lz_u * ( Den(t-1)*Cap(t-1)*x1*y2 + ...
                                                    Den(t)*Cap(t)*x2*y2);   
                        value(pointer, 3) = h*lx*lz+h_d*lx*ly;
                        pointer=pointer+1;                      

                end
            end
        end
    end
    fprintf('package bot done: %f\n', toc(t_check));
    fprintf('%i of %i layers are analyzed\n', checker, Layer.Nt+1)
    toc(t_start);
    fprintf('\n');
    row = row(1:pointer-1,1);
    column = column(1:pointer-1,1);
    value = value(1:pointer-1,:);
    Y = sparse(row, column, value(:,1), var, var);
    C = sparse(row, column, value(:,2), var, var);
    H = sparse(row, column, value(:,3), var, var);
end

