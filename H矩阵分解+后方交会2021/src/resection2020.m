clc;
clear;

format long g
x0=1487.240000;y0=1973.510000;
f=(3029.350000+3032.070000)/2;
k1=0.0488249;k2=-0.168801;k3=0.194688;p1=-0.000276;p2=0.044539;

m=[1357 811
   1366 2359
   2201 2426
   2210 635];

M=[0.39 1.213 1.960
0.39 1.213 0.965
0.948 1.213 0.965
0.948 1.213 1.960];


[row,col]=size(M);

%step1 对二维像点坐标去畸变并标准化像素
un_m=zeros(row,2);
for i=1:row
%     [un_x,un_y]=removeDistortion(m(i,1),m(i,2),x0,y0,f,k1,k2,k3,p1,p2);
    un_m(i,1)=(m(i,1)-x0)/f;
    un_m(i,2)=(m(i,2)-y0)/f;
end

%step2 对三维空间点坐标进行重心化并对其进行svd分解
Mc=sum(M)/row;
MM=zeros(row,3);
for i=1:row
    MM(i,:)=M(i,:)-Mc;
end
mm=MM'*MM;
[U,W,V]=svd(mm);

%step3 将svd分解的结果作为子旋转矩阵，并将三维空间坐标旋转平移至xoy平面上
R_temp=V';
T_temp=-R_temp*Mc';
Mxy=zeros(row,2);
for i=1:row
    XYm=zeros(row,2);
    XYm(i,1)=R_temp(1,:)*M(i,:)'+T_temp(1);
    XYm(i,2)=R_temp(2,:)*M(i,:)'+T_temp(2);
    Mxy(i,:)=[XYm(i,1),XYm(i,2)];
end

%计算H-矩阵
H = calH(Mxy,un_m)
% H=[-0.566, 0.0171,0.1772
%     -0.03102,-0.6684,-0.1887
%     0.3443,-0.03058,1];
h1=H(:,1);h2=H(:,2);t=H(:,3);

%分解H-矩阵获得相机位姿
h1_norm=sqrt(h1(1)*h1(1)+h1(2)*h1(2)+h1(3)*h1(3));
h2_norm=sqrt(h2(1)*h2(1)+h2(2)*h2(2)+h2(3)*h2(3));
r1=h1/h1_norm;
r2=h2/h2_norm;
t=t/((h1_norm+h2_norm)/2);
r3=cross(r1,r2);
r3=r3/norm(r3);
R=[r1,r2,r3]
T=R * T_temp + t
R_last=R*R_temp
T_last=-inv(R_last)*T %计算线元素（摄影中心坐标）
 
 [r1,r2,r3] = dcm2angle(R_last,'XYZ');

sum=0;
diffEO = 1;
EO=[T_last(1),T_last(2),T_last(3),r1,r2,r3]';
xi=zeros(row,3);
XYZm=zeros(row,3);
xy_temp=zeros(row,2);
A=zeros(row,6);
L=zeros(row,1);
%% omega phi kappa
while diffEO > 1e-5 && sum < 500

    omega = EO(4);
    phi = EO(5);
    kappa = EO(6);
    C1 = EO(1:3);
    
    a1=cos(phi)*cos(kappa);
    a2=-cos(phi)*sin(kappa);
    a3=sin(phi);
    b1=sin(omega)*sin(phi)*cos(kappa)+cos(omega)*sin(kappa);
    b2=-sin(omega)*sin(phi)*sin(kappa)+cos(omega)*cos(kappa);
    b3=-sin(omega)*cos(phi);
    c1=-cos(omega)*sin(phi)*cos(kappa)+sin(omega)*sin(kappa);
    c2=cos(omega)*sin(phi)*sin(kappa)+sin(omega)*cos(kappa);
    c3=cos(omega)*cos(phi);
    R=[a1 b1 c1;
       a2 b2 c2;
       a3 b3 c3];
   
   cosw = cos(omega);
   sinw = sin(omega);  
   cosp = cos(phi);
   sinp = sin(phi);
   cosk = cos(kappa);
   sink = sin(kappa);

    for i=1:size(m,1)
        xi(i,:) = M(i,1:3)'- C1; %XA-XS
        XYZm(i,:) = R * xi(i,:)';
        
        xy_temp(i,1) = f * XYZm(i,1) / XYZm(i,3) + x0;
        xy_temp(i,2) = f * XYZm(i,2)/XYZm(i,3) + y0;
        xyp = xy_temp(i,:);
        
        X_ = XYZm(i,1);
        Y_ = XYZm(i,2);
        Z_ = XYZm(i,3);
        
        x_x0 = xyp(1) - x0;
        y_y0 = xyp(2) - y0;
        
        X_XS = xi(i,1);
        Y_YS = xi(i,2);
        Z_ZS = xi(i,3);
        
        dx_XS = 1 / Z_ * (-f * a1 + a3 * x_x0);
        dx_YS = 1 / Z_ * (-f * b1 + b3 * x_x0);
        dx_ZS = 1 / Z_ * (-f * c1 + c3 * x_x0);
        dx_omega = f / Z_ * (b1 * Z_ZS- c1 * Y_YS) - x_x0 / Z_ * (b3 * Z_ZS - c3 * Y_YS);       
        temp_one = f / Z_ * (-sinp*cosk * X_XS + sinw*cosp*cosk * Y_YS - cosw*cosp*cosk * Z_ZS);
        temp_two = -x_x0 / Z_ * (cosp * X_XS + sinw*sinp * Y_YS - cosw*sinp * Z_ZS);
        dx_phi =  temp_one + temp_two;               
        dx_kappa = f / Z_ * Y_;   
%         dx_kappa = y_y0 ;
        
        dy_XS = 1 / Z_ * (-f * a2 + a3 * y_y0);
        dy_YS = 1 / Z_ * (-f * b2 + b3 * y_y0);
        dy_ZS = 1 / Z_ * (-f * c2 + c3 * y_y0);
        dy_omega = f / Z_ * (b2 * Z_ZS- c2 * Y_YS) - y_y0 / Z_ * (b3 * Z_ZS - c3 * Y_YS);        
        temp_one = f / Z_ * (sink*sinp * X_XS - sinw*sink*cosp * Y_YS + cosw*sink*cosp * Z_ZS);
        temp_two = -y_y0 / Z_ * (cosp * X_XS + sinw*sinp * Y_YS - cosw*sinp * Z_ZS);
        dy_phi =  temp_one + temp_two;
        dy_kappa = -f / Z_ * X_; 
%         dy_kappa = -x_x0 ;
        
        A(2*i-1,:)=[dx_XS dx_YS dx_ZS dx_omega dx_phi dx_kappa];                             
        A(2*i  ,:)=[dy_XS dy_YS dy_ZS dy_omega dy_phi dy_kappa]; 
        L(2*i-1)=m(i,1)-xyp(1);
        L(2*i)=m(i,2)-xyp(2);
    end
 
    dEO = inv(A'*A) * (A'*L);
    diffEO = max(dEO);
    EO = EO + dEO;
    sum = sum + 1;
end
C1
sum
R
