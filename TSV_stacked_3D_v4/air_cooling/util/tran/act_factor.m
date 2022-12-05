%rng(11223344);

pulse_width_l = 0.1;
pulse_width_u = 0.1;
len = 100;
max_act = 0.8;
min_act = 0;
array = zeros(5000, 1);
act = zeros(5000,1);

Tcurrent = 0;
l = 1;
while(len - Tcurrent > pulse_width_u)
    array(l) = Tcurrent;
    act(l) = min_act + (max_act - min_act)*rand(1);
    temp = pulse_width_l + (pulse_width_u - pulse_width_l)*rand(1);    
    l = l + 1;
    Tcurrent = Tcurrent + temp;
end
array(l) = len;
act(l) = min_act + (max_act - min_act)*rand(1);

array = array(1:l);
act = act(1:l);

%stairs(array, act);

fid=fopen('time.txt','w+');
fprintf(fid,'%0.1f\r\n',array);
fclose(fid);

fid=fopen('act.txt','w+');
fprintf(fid,'%0.2f\r\n',act);
fclose(fid);

array = load('time.txt');
act = load('act.txt');
stairs(array, act);

