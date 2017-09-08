clear all
close all

load('oneVesselFenics.mat');
saveasfilename = 'oneVesselGrid_d10';
d = 10;

% Stop doing things
hole_coor = uint8(hole_coor./d)*d;
Hx = double(corners(1,1)):d:double(corners(2,1));
Hy = double(corners(1,2)):d:double(corners(2,2));
[X, Y] = meshgrid(Hx, Hy);
r = cell(1,2);
r1 = sqrt((X-double(hole_coor(1,1))).^2 + (Y-double(hole_coor(1,2))).^2);
r2 = sqrt((X-double(hole_coor(2,1))).^2 + (Y-double(hole_coor(2,2))).^2);
r{1} = r1;
r{2} = r2;

Nx = length(Hx);
Ny = length(Hy);
Nxy = length(mesh_coor);
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
                xs = mesh_coor(k,1);
                ys = mesh_coor(k,2);
                xixs = abs(xi-xs);
                yjys = abs(yj-ys);
                totdifference = xixs + yjys;
                if totdifference < mytotdifference
                    myindex = k;
                    mytotdifference = totdifference;                
                end            
            end
            P(i,j) = p_array(myindex);
        end
    end
end

% Dette tar vekk for mye. Du må trekke fra gjennomsnittet av alle SE.
% Hvordan fant vi 3? Den egentlig verdien er 3.2754. Se i notater.
% mydiff = p_noisy - Psm.
% Fjern største verdi n_åre ganger. Dette er en outlier.
% std(mydiff(:)) = 3.2754
%P_noisy = normrnd(P, 3.2754-mean(SE(:)));
%rng(seed2);
%P_noisy = normrnd(P, SE);

rng(1);
P_noisy = normrnd(P, 3);

figure(1)
imagesc(Hx, Hy, P, [0, max(P(:))])
title('\textbf{Ground truth model data}', 'Interpreter', 'latex')
xlabel(['$x\, [',num2str(d),' \mu m]$'], 'Interpreter', 'latex');
ylabel(['$y\, [',num2str(d),' \mu m]$'], 'Interpreter', 'latex');   
colormap(makeColorMap([1,1,1], [1,0,0], 1000));
h = colorbar;
xlabel(h,'$\mathrm{pO_2}$ [mmHg]', 'Interpreter', 'latex')
set(gca, 'fontsize', 16);   

figure(2)
imagesc(Hx, Hy, P_noisy, [0, max(P(:))])
title('\textbf{Ground truth model data w/ noise}', 'Interpreter', 'latex')
xlabel(['$x\, [',num2str(d),' \mu m]$'], 'Interpreter', 'latex');
ylabel(['$y\, [',num2str(d),' \mu m]$'], 'Interpreter', 'latex');   
colormap(makeColorMap([1,1,1], [1,0,0], 1000));
h = colorbar;
xlabel(h,'$\mathrm{pO_2}$ [mmHg]', 'Interpreter', 'latex')
set(gca, 'fontsize', 16); 

save(saveasfilename, 'P', 'P_noisy', 'r', 'd', 'Hx', 'Hy', 'M_true', 'p_ves', 'r_ves', 'corners', 'hole_coor', 'Nx', 'Ny', 'Nxy')
