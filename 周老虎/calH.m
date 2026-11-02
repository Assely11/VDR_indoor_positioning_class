function [ H ] = calH( A,B )
%CALH 此处显示有关此函数的摘要
%   此处显示详细说明

% % AX = B
% A=[];
% B=[];
 [m,~]=size(A);
 M=zeros(m,8);
 N=zeros(m,1);
 for i=1:m
     u1 = A(i,1);    v1 = A(i,2);
     u2 = B(i,1);    v2 = B(i,2);
     M(2*i-1,:)=[  u1    v1    1.0    0.0   0.0    0.0   -u1*u2    -v1*u2];
     M(2*i,:)=[  0.0   0.0   0.0     u1    v1    1.0   -u1*v2   -v1*v2 ];
     N(2*i-1,:)=(u2);
     N(2*i,:)=(v2);
 end

 X=inv((M'*M))*M'*N;
 H = [X(1) X(2) X(3) 
     X(4) X(5) X(6) 
     X(7) X(8) 1.0]; 
end

