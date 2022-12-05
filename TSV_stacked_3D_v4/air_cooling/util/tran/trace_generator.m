%%%%%chip trace generator%%%%%%%%%%%%%%%%
FLAG_GRADUAL = 0;
power = 74.49;
map = [ 0.2e-3 0.2e-3 1.0e-3 4.6e-3 2.3000;
        0.2e-3 5.0e-3 1.0e-3 3.6e-3 2.2690;
        1.4e-3 0.2e-3 7.3e-3 2.6e-3 14.8044;
        1.4e-3 3.0e-3 1.4e-3 5.6e-3 7.5264;
        3.0e-3 3.0e-3 1.3e-3 5.6e-3 8.3720;
        4.5e-3 3.0e-3 1.0e-3 3.5e-3 1.9600;     
        4.5e-3 6.7e-3 1.0e-3 1.9e-3 1.3680;
        5.7e-3 3.0e-3 1.4e-3 5.6e-3 6.9776;
        7.3e-3 3.0e-3 1.4e-3 5.6e-3 10.5840;
        8.9e-3 0.2e-3 1.0e-3 4.2e-3 2.4360;
        8.9e-3 4.6e-3 1.0e-3 4.0e-3 3.0800;
        0.2e-3 8.8e-3 9.7e-3 1.0e-3 6.4020];
%trace format:
%first line: time
%second line: how many blocks #B
%third line: total power
%the following #B lines list the power of each block
time = 1;
ratio = 1;
array = load('time.txt');
act = load('act.txt');
act(1)=0.48;
figure();
stairs(array*ratio, act);
len = length(array);

fid=fopen('Chip1.trace','w+');
Trise = ratio*0.1;
Tfall = ratio*0.05;
Tstep = ratio*0.01;

Nrise = Trise/Tstep+1;
Nfall = Tfall/Tstep+1;

l = 1;
t = 0;
x = [t*ratio 12 power*act(l)];
fprintf(fid,'%g\r\n', x);

x = map(:,5)*act(1);
fprintf(fid,'%g\r\n', x);

for l = 1:1:len-1
    if FLAG_GRADUAL == 1
        if act(l+1) > act(l)
            t_next = array(l+1)*ratio - Trise;
            temp = linspace(act(l), act(l+1), Nrise);
            N = Nrise;
        elseif act(l+1) < act(l)
            t_next = array(l+1)*ratio - Tfall;
            temp = linspace(act(l), act(l+1), Nfall);
            N = Nfall;
        else
            N = 0;
        end
        for j = 1:N
            x = [t_next+Tstep*(j-1), 12, power*temp(j)];
            fprintf(fid,'%g\r\n', x);
            x = map(:,5)*temp(j);
            fprintf(fid,'%g\r\n', x);
        end
    else
            x=[Tstep*l, 12, power*act(l+1)];
            fprintf(fid,'%g\r\n', x);
            x = map(:,5)*act(l+1);
            fprintf(fid,'%g\r\n', x);
    end
end

fclose(fid);
    
    