clc;
clear;

XYZw = [1 0 0
        0 1 0
        0 0 1.0];

tpl  = mean(XYZw)';
XYZc = XYZw(:,:) - tpl';
[U,W,V] = svd(XYZc'*XYZc);
Rpl = V;

Tpl = Rpl;
Tpl(:,4) = tpl;
Tpl = [Tpl ; 0 0 0 1];

XYZw(:,4) = 1.0;
XYZpl = inv(Tpl) * XYZw';
XYZpl = XYZpl';
XYZpl = XYZpl(:,1:3);

%%

XYZw(4,:) = XYZw(1,:);
XYZpl(4,:) = XYZpl(1,:);

plot3(XYZw(:,1),XYZw(:,2),XYZw(:,3),'r');
hold on
plot3(XYZpl(:,1),XYZpl(:,2),XYZpl(:,3),'g');

