function [ undistor_x,undistor_y ] = removeDistortion( distor_x,distor_y,x0,y0,f,k1,k2,k3,p1,p2 )
%REMOVEDISTORTION 此处显示有关此函数的摘要
%   此处显示详细说明
   undistor_x=distor_x;
   undistor_y=distor_y;
   error = 1;
   
   while error >1e-10
   xd=(undistor_x-x0)/f;
   yd=(undistor_y-y0)/f;   
   r2=xd*xd+yd*yd;
   r4=r2*r2;
   r6=r4*r2;
   x2=xd*xd;
   y2=yd*yd;
   xy=xd*yd;   
   dx=xd*(k1*r2+k2*r4+k3*r6)+2*p1*xy+p2*(r2+2*x2);
   dy=yd*(k1*r2+k2*r4+k3*r6)+p1*(r2+2*y2)+2*p2*xy;
      distor_x_=xd+dx;
   distor_y_=yd+dy;   
   distor_x_=distor_x_*f+x0;
   distor_y_=distor_y_*f+y0;   
   error = abs(distor_x-distor_x_)+abs(distor_y_-distor_y_);   
   undistor_x=distor_x-dx*f;
   undistor_y=distor_y-dy*f;
   
   end

end

