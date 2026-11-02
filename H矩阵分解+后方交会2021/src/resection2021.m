clc;
clear;

format long g
dataPathName = '.\..\data_三层走廊\mate20\位置4';
innerPara = load([dataPathName, '\innerPara.txt']);
XYZw = load([dataPathName, '\objPoints.txt']);
uv_pix = load([dataPathName, '\imgPoints.txt']);

f = innerPara(1);
u0 = innerPara(2);
v0 = innerPara(3);
k1  = innerPara(4); k2  = innerPara(5); k3  = innerPara(6);
p1  = innerPara(7); p2  = innerPara(8); 


K = [f 0  u0
    0  f  v0
    0  0 1];


pointNum = size(uv_pix,1);

%step1 对二维像点坐标去畸变并标准化像素
for i = 1:pointNum
 % [uv_pix(i,1),uv_pix(i,2)] = removeDistortion( uv_pix(i,1),uv_pix(i,2),u0,v0,f,k1,k2,k3,p1,p2 );
end
uv_pix(:,3) = 1;
xy_np = inv(K) * uv_pix';


%step2 对三维空间点坐标进行重心化并对其进行svd分解
%      将svd分解的结果作为子旋转矩阵，并将三维空间坐标旋转平移至xoy平面上
tpl  = mean(XYZw)';
XYZc = XYZw(:,:) - tpl';
[U,W,V] = svd(XYZc'*XYZc);
Rpl = V;

Tpl = Rpl;
Tpl(:,4) = tpl;
Tpl = [Tpl ; 0 0 0 1];

XYZw(:,4) = 1.0;
XYZpl = inv(Tpl) * XYZw';

%step3 计算H-矩阵,分解H-矩阵获得相机位姿
XYZpl = XYZpl';
%XYZpl = XYZpl(:,1:2);
xy_np = xy_np';
%xy_np = xy_np(:,1:2);
H = calH(XYZpl,xy_np);


r1 = H(:,1)/norm(H(:,1));
r2 = H(:,2)/norm(H(:,2));
r3 = cross(r1,r2);
r2 = cross(r3,r1);
Rh = [r1,r2,r3];
R = Rh*Rpl';

lamda = (norm(H(:,1))+norm(H(:,2))) / 2;
th = H(:,3) / lamda;
C = -R'*(th - R * tpl);

[omega phi kappa] = dcm2angle(R,'XYZ');
OPK = [omega phi kappa];
fprintf('C 初始值： %f\t%f\t%f\n', C(1),C(2),C(3));
fprintf('OPK 初始值： %f\t%f\t%f\n', omega*180/pi, phi*180/pi, kappa*180/pi);

%R1 = angle2dcm(omega,phi, kappa,'XYZ')
%R2 = opk2R(omega,phi, kappa)

%step4 优化EO
[betterR  betterC] = resection_iteration(K, R, C,  XYZw, uv_pix);
[omega phi kappa] = dcm2angle(betterR ,'XYZ');

fprintf('\n');
fprintf('C 优化值： %f\t%f\t%f\n', betterC(1),betterC(2),betterC(3));
fprintf('OPK 优化值： %f\t%f\t%f\n', omega*180/pi, phi*180/pi, kappa*180/pi);
fprintf('done!\n');

%
% C 初始值： 1.414088	-0.483548	1.289718
% OPK 初始值： -92.555121	-28.707038	-0.986163
% 
% C 优化值： 1.447536	-0.454899	1.277983
% OPK 优化值： -92.161524	-29.990539	-0.853488
%