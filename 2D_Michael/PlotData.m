umax_1=1;
umax_2=100;
[u1_grid,u2_grid] = meshgrid(0:0.1:umax_1,0:10:umax_2);
SSAData=cell(size(u1_grid,1),size(u1_grid,2));
xyz = [];
for i=1:size(u1_grid,1)
    for j=1:size(u2_grid,2)
        s = dlmread(strcat('2D_Michael/Data/SSAData_',num2str(i),'_',num2str(j)),'\t');
        xyz = [xyz; [u1_grid(i,j)*ones(1,20)' u2_grid(i,j)*ones(1,20)' s(:,1)]];
    end
end
hold on
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'.b')
hold off