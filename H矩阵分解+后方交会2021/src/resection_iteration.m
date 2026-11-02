function [betterR, betterC] = resection_iteration(K, R, C,  XYZw, uv_pix)

pointNum = size(XYZw,1);
f = K(1,1);
x0 = K(1,3);
y0 = K(2,3);
[omega phi kappa] = dcm2angle(R,'XYZ');



sum=0;
diffEO = 1;
EO=[C(1),C(2),C(3),omega,phi,kappa]';
xi=zeros(pointNum,3);
XYZm=zeros(pointNum,3);
xy_temp=zeros(pointNum,2);
A=zeros(pointNum,6);
L=zeros(pointNum,1);
%% omega phi kappa
while diffEO > 1e-5 && sum < 100

    omega = EO(4);
    phi = EO(5);
    kappa = EO(6);
    C1 = EO(1:3);
        
    %R = opk2R(omega,phi, kappa);
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

    for i=1:pointNum
        xi(i,:) = XYZw(i,1:3)'- C1; %XA-XS
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
        L(2*i-1)=uv_pix(i,1)-xyp(1);
        L(2*i)=uv_pix(i,2)-xyp(2);
    end
 
    dEO = inv(A'*A) * (A'*L);
    diffEO = max(dEO);
    EO = EO + dEO;
    sum = sum + 1;
end
    betterR =  opk2R(EO(4), EO(5), EO(6));
    betterC = EO(1:3);
end

