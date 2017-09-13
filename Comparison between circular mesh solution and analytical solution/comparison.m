clear all
close all

% ******************************************
% P_num and P_anal
% ******************************************

load('circleMesh_res200_d1.mat')
P_anal = 80 + 0.25*M_true*(r.^2 - 6^2) - 0.5*M_true*200^2*log(r./6);
P_anal(r < 6) = 80;
difference = P - P_anal;

figure(1)
imagesc(Hx, Hy, P, [0, max(P(:))]);
title('$\mathrm{P_{num}}$', 'Interpreter', 'latex');
xlabel('$x\, [\mu m]$', 'Interpreter', 'latex');
ylabel('$y\, [\mu m]$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
colormap(makeColorMap([1,1,1], [1,0,0], 1000));
h = colorbar; axis xy;
xlabel(h,'$\mathrm{P_{num}}$', 'Interpreter', 'latex')

figure(2)
imagesc(Hx, Hy, P_anal, [0, max(P_anal(:))]);
title('$\mathrm{P_{anal}}$', 'Interpreter', 'latex');
xlabel('$x\, [\mu m]$', 'Interpreter', 'latex');
ylabel('$y\, [\mu m]$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
colormap(makeColorMap([1,1,1], [1,0,0], 1000));
h = colorbar; axis xy;
xlabel(h,'$\mathrm{P_{anal}}$', 'Interpreter', 'latex')

figure(3)
imagesc(Hx, Hy, difference);
title('$\mathrm{P_{num}} - \mathrm{P_{anal}}$', 'Interpreter', 'latex');
xlabel('$x\, [\mu m]$', 'Interpreter', 'latex');
ylabel('$y\, [\mu m]$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
colormap(NegativeEnhancingColormap(1000, [min(difference(:)) max(difference(:))], [0 0 1], [1 0 0], 1));
h = colorbar; axis xy;

% ******************************************
% Second derivatives.
% Ground truth second derivative is M_true.
% ******************************************

del2P = 4*del2(P, double(d));
del2P_anal = 4*del2(P_anal, double(d));
epsilon = (del2P - M_true) ./ M_true;
epsilon_anal = (del2P_anal - M_true) ./ M_true;

figure(4);        
del2P(r < 20) = 0;             
imagesc(Hx, Hy, del2P);
title('\texttt{del2P(P)}', 'Interpreter', 'latex');
xlabel('$x\, [\mu m]$', 'Interpreter', 'latex');
ylabel('$y\, [\mu m]$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
colormap(NegativeEnhancingColormap(1000, [min(del2P(:)) max(del2P(:))], [0 0 1], [1 0 0], 1));
h = colorbar; axis xy;
xlabel(h, '\texttt{del2P(P)}', 'Interpreter', 'latex');

figure(5);        
del2P_anal(r < 20) = 0;             
imagesc(Hx, Hy, del2P_anal);
title('\texttt{del2P(P\_anal)}', 'Interpreter', 'latex');
xlabel('$x\, [\mu m]$', 'Interpreter', 'latex');
ylabel('$y\, [\mu m]$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
colormap(makeColorMap([1,1,1], [1,0,0], 1000));
h = colorbar; axis xy;
xlabel(h, '\texttt{del2P(P\_anal)}', 'Interpreter', 'latex');

figure(6);
epsilon(r < 20) = 0;             
imagesc(Hx, Hy, epsilon);   
title(['$(\texttt{del2P(P)}-\mathrm{M_{true}})/\mathrm{M_{true}}$'], 'Interpreter', 'latex'); 
xlabel('$x\, [\mu m]$', 'Interpreter', 'latex');
ylabel('$y\, [\mu m]$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
colormap(NegativeEnhancingColormap(1000, [min(epsilon(:)) max(epsilon(:))], [0 0 1], [1 0 0], 1));
h = colorbar; axis xy;
xlabel(h, 'epsilon', 'Interpreter', 'latex');

figure(7);
epsilon_anal(r < 20) = 0;             
imagesc(Hx, Hy, epsilon_anal);   
title(['$(\texttt{del2P(P\_anal)}-\mathrm{M_{true}})/\mathrm{M_{true}}$'], 'Interpreter', 'latex'); 
xlabel('$x\, [\mu m]$', 'Interpreter', 'latex');
ylabel('$y\, [\mu m]$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
colormap(NegativeEnhancingColormap(1000, [min(epsilon_anal(:)) max(epsilon_anal(:))], [0 0 1], [1 0 0], 1));
h = colorbar; axis xy;
xlabel(h, '$\mathrm{epsilon_{anal}}$', 'Interpreter', 'latex');


% ******************************************
% Resolution
% ******************************************

resolution = zeros(1,7);

minDifference_d1 = zeros(1,7);
maxDifference_d1 = zeros(1,7);
meanDifference_d1 = zeros(1,7);

minDifference_d10 = zeros(1,7);
maxDifference_d10 = zeros(1,7);
meanDifference_d10 = zeros(1,7);

for i = 2:9
    res = 100*i;
    filename = ['circleMesh_res', num2str(res), '_d1.mat'];
    load(filename)
    P_anal = 80 + 0.25*M_true*(r.^2 - 6^2) - 0.5*M_true*200^2*log(r./6);
    P_anal(r < 6) = 80;
    difference = P - P_anal;
    resolution(i-1) = res;
    minDifference_d1(i-1) = min(abs(difference(:)));
    maxDifference_d1(i-1) = max(abs(difference(:)));
    meanDifference_d1(i-1) = mean(abs(difference(:)));
    
    filename = ['circleMesh_res', num2str(res), '_d10.mat'];
    load(filename)
    P_anal = 80 + 0.25*M_true*(r.^2 - 6^2) - 0.5*M_true*200^2*log(r./6);
    P_anal(r < 6) = 80;
    difference = P - P_anal;    
    minDifference_d10(i-1) = min(abs(difference(:)));
    maxDifference_d10(i-1) = max(abs(difference(:)));
    meanDifference_d10(i-1) = mean(abs(difference(:)));
end

figure(8)
plot(resolution, meanDifference_d1, 'r-s')
hold on
plot(resolution, meanDifference_d10, 'b-o')
title('\texttt{mean(abs(P\_num - P\_anal))}', 'Interpreter', 'latex')
xlabel('Resolution')
ylabel('Mean')
legend('d = 1', 'd = 10')
set(gca, 'fontsize', 16);

figure(9)
plot(resolution, minDifference_d1, 'r-s')
hold on
plot(resolution, minDifference_d10, 'b-o')
title('\texttt{min(abs(P\_num - P\_anal))}', 'Interpreter', 'latex')
xlabel('Resolution')
ylabel('Min')
legend('d = 1', 'd = 10')
set(gca, 'fontsize', 16);

figure(10)
plot(resolution, maxDifference_d1, 'r-s')
hold on
plot(resolution, maxDifference_d10, 'b-o')
title('\texttt{max(abs(P\_num - P\_anal))}', 'Interpreter', 'latex')
xlabel('Resolution')
ylabel('Max')
legend('d = 1', 'd = 10')
set(gca, 'fontsize', 16);