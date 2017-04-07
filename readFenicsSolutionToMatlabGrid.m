clear all
close all

load('po2grid.mat');
d = 5;
filename = 'testgrid';

% Stop doing things
x0 = corners(1);
y0 = corners(2);
x1 = corners(3);
y1 = corners(4);

Hx = x0:d:x1;
Hy = y0:d:y1;
[X, Y] = meshgrid(Hx, Hy);
r = cell(1,2);
r1 = sqrt((X-double(center1(1))).^2 + (Y-double(center1(2))).^2);
r2 = sqrt((X-double(center2(1))).^2 + (Y-double(center2(2))).^2);
r{1} = r1;
r{2} = r2;

Nx = length(Hx);
Ny = length(Hy);
Nxy = length(xyvec);
P = NaN(Ny, Nx);

for i = 1:Ny
    for j = 1:Nx
        xi = Hx(j);
        yj = Hy(i);
        if (r1(i,j) < r_ves)
            P(i,j) = p_ves(1);
        elseif (r2(i,j) < r_ves)
            P(i,j) = p_ves(2);
        else
            myindex = 0;
            mytotdifference = 10000;
            for k = 1:Nxy
                xs = xyvec(k,1);
                ys = xyvec(k,2);
                xixs = abs(xi-xs);
                yjys = abs(yj-ys);
                totdifference = xixs + yjys;
                if totdifference < mytotdifference
                    myindex = k;
                    mytotdifference = totdifference;                
                end
            P(i,j) = cvec(myindex);
            end
        end
    end
end

imagesc(Hx, Hy, P, [0, max(P(:))])
title('\textbf{Ground truth model data}', 'Interpreter', 'latex')
xlabel(['$x\, [',num2str(d),' \mu m]$'], 'Interpreter', 'latex');
ylabel(['$y\, [',num2str(d),' \mu m]$'], 'Interpreter', 'latex');   
colormap(makeColorMap([1,1,1], [1,0,0], 1000));
h = colorbar;
xlabel(h,'$\mathrm{pO_2}$ [mmHg]', 'Interpreter', 'latex')
set(gca, 'fontsize', 16);   

save(filename, 'd', 'Hx', 'Hy', 'P', 'r', 'M_true', 'corners')