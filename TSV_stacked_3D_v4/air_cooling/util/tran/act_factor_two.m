rng(11223344);

pulse_width_l = 0.1;
pulse_width_u = 0.1;
len = 100;
max_act = 0.78;
min_act = 0.01;
array = zeros(1000, 1);
act = zeros(1000,1);

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

array_short = array(1:10:end);
act_short = zeros(size(array_short));
for i = 1:length(act_short)
    END = min(10+(i-1)*10, size(act,1));
    act_short(i) = mean(act(1+(i-1)*10:END));
end

%stairs(array, act);

fid=fopen('time.txt','w+');
fprintf(fid,'%0.1f\r\n',array);
fclose(fid);

fid=fopen('act.txt','w+');
fprintf(fid,'%0.2f\r\n',act);
fclose(fid);

fid=fopen('time_short.txt','w+');
fprintf(fid,'%0.1f\r\n',array_short);
fclose(fid);

fid=fopen('act_short.txt','w+');
fprintf(fid,'%0.2f\r\n',act_short);
fclose(fid);

array = load('time.txt');
act = load('act.txt');
h2 = stairs(array, act, 'color', 'b', 'LineWidth', 2);
hold on;

% Trange = [0 1];
% set(gca, 'ylim', Trange, 'YTick', linspace(Trange(1), Trange(2), 10));

array_short = load('time_short.txt');
act_short = load('act_short.txt');
h1 = stairs(array_short, act_short, 'color', 'r', 'LineWidth', 3);

set(gca, 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Time(s)');
ylabel('Activity');

% LE1 = 'Coarse Sampling';
% LE2 = 'Fine Sampling';
% 
% lgnd = legend([h1, h2], LE1, LE2);
% set(lgnd,'color','none');
