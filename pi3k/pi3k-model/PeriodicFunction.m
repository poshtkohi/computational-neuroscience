maxNumCompThreads(1);
%parpool(2);
%delete(gcp('nocreate'));
%----- Prepare the Environment -------%
close all; clc; clear; clearvars;
format long g;
TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
%----- simulation -------%
t_max = 10; Z% integer
n = 1000;
x = zeros(n, 1);
t = linspace(0, t_max, n);
period = 1;
for i=1:1:length(x)
    x(i) = periodic_function(t(i), period);
end

x(1) = 0;

fig1 = figure('Name', 'Periodic Function');
fig1.WindowStyle = 'docked';

hold on;
plot(t, x, '-ko', 'MarkerIndices', 1:2:length(x));
hold off;
%legend(sprintf('[Ca^2+_i] at ADP=%gmM', ATP(1)), sprintf('[Ca^2+_i] at ADP=%gmM', ATP(2)), sprintf('[Ca^2+_i] at ADP=%gmM', ATP(3)), sprintf('[Ca^2+_i] at ADP=%gmM', ATP(4)), sprintf('[Ca^2+_i] at ADP=%gmM', ATP(5)), 'FontWeight', 'bold');
xlabel('Time (s)');
%ylabel('[Ca^2+_i] (nM)');
%xlim([0 t_max]);
%------Functions -------%
function [val] = periodic_function(t, period)
   num = ceil(t);
   if mod(num, 2) == 0
       val = 1;
   else
       val = 0;
   end
   t
end
%--------------------%