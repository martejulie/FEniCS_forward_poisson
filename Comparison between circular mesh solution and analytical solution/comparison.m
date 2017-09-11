
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
    po2_anal = 80 + 0.25*M_true*(r.^2-6^2) - 0.5*M_true*200^2*log(r./6);
    po2_anal(r<6)=80;
    difference = P-po2_anal;
    resolution(i-1) = res;
    minDifference_d1(i-1) = min(abs(difference(:)));
    maxDifference_d1(i-1) = max(abs(difference(:)));
    meanDifference_d1(i-1) = mean(abs(difference(:)));
    
    filename = ['circleMesh_res', num2str(res), '_d10.mat'];
    load(filename)
    po2_anal = 80 + 0.25*M_true*(r.^2-6^2) - 0.5*M_true*200^2*log(r./6);
    po2_anal(r<6)=80;
    difference = P-po2_anal;    
    minDifference_d10(i-1) = min(abs(difference(:)));
    maxDifference_d10(i-1) = max(abs(difference(:)));
    meanDifference_d10(i-1) = mean(abs(difference(:)));
end

figure(1)
plot(resolution, meanDifference_d1, 'r-s')
hold on
plot(resolution, meanDifference_d10, 'b-o')
title('Mean difference')
xlabel('Resolution')
ylabel('Mean difference')

figure(2)
plot(resolution, maxDifference_d1, 'r-s')
hold on
plot(resolution, maxDifference_d10, 'b-o')
title('Max difference')
xlabel('Resolution')
ylabel('Max difference')