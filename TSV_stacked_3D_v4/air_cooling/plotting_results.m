clear;
clc;
close all;

ratio = 0.1;
array = load('time.txt');
act = load('act.txt');
[xst, yst] = stairs(array*ratio, act);

time = load('Tfinal.txt');
Tmax_pro = load('./result/tran/chip1_die1_max.txt');
Tmax_fpga = load('./result/tran/chip2_die1_max.txt');
Tmax_mem = load('./result/tran/chip3_die5_max.txt');

% Tmax_pro(end+1) = Tmax_pro(end);
% Tmax_fpga(end+1) = Tmax_fpga(end);
% Tmax_mem(end+1) = Tmax_mem(end);

 Trange = [64 88];
 Prange = [0 70];

figure(3);
[ax, l1, l2] = plotyy(time,Tmax_pro, xst, yst*74.49);
hold(ax(1));

t_line = 3;
p_line = 2;

set(l1, 'LineStyle', '--', 'linewidth', t_line, 'color', 'r')
set(l2, 'linewidth', p_line, 'color', 'k')
set(ax(2),'Ycolor','k')
set(ax(1),'Ycolor','k')
set(ax(1), 'ylim', Trange, 'YTick', linspace(Trange(1), Trange(2), 7));
set(ax(2), 'ylim', Prange, 'ytick', linspace(Prange(1), Prange(2), 6));
%grid on;
set(ax(1), 'FontSize', 16, 'FontWeight', 'bold');
set(ax(2), 'FontSize', 16, 'FontWeight', 'bold');
xlabel(ax(1), 'Time(s)');
ylabel(ax(1), 'Maximum T_{junc} (', char(0176), 'C)');
ylabel(ax(2), 'Processor Power (W)');
h1 = plot(ax(1), time, Tmax_fpga);
set(h1, 'LineStyle', '--','linewidth', t_line, 'color', 'k');
h2 = plot(ax(1), time, Tmax_mem);
set(h2, 'LineStyle', '--','linewidth', t_line, 'color', 'b');

LE1 = 'CPU';
LE2 = 'FPGA';
LE3 = 'DRAM';
LE4 = 'Processor Power';

lgnd = legend([l1, h1, h2, l2], LE1, LE2, LE3, LE4);

set(lgnd,'color','none');