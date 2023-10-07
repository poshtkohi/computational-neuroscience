maxNumCompThreads(1);
%parpool(2);
%delete(gcp('nocreate'));
%----- Prepare the Environment -------%
close all; clc; clear; clearvars;
format long g;
TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
%----- simulation -------%
global t1 t2 g enabled;
g = 0.5;
enabled = 1;
t1 = 2;
t2 = 4;
t_max = 20;
delay = 2;
width = t2 - t1;
amplitude = 50;

t = linspace(0, t_max, 1000);
x = zeros(length(t), 1);
G = ones(length(t), 1);
 
for i=1:1:length(t)
    [x(i), G(i)] = periodic_function(t(i), delay, width, amplitude);
end

fig1 = figure('Name', 'Periodic Function');
fig1.WindowStyle = 'docked';

hold on;
%x = normalize(x, 'range');
%G = normalize(G, 'range');
plot(t, x, '-ko', 'MarkerIndices', 1:50:length(x));
plot(t, G, 'b--*', 'MarkerIndices', 1:50:length(G));
hold off;
legend('ADP', 'g', 'FontWeight', 'bold');
xlabel('Time (s)');
%--------------------%