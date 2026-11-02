clc;
clear;
format long g
img = imread('ºŸ¿œª¢.jpg');
imgPoint = [284 115.5
            274 166 
            179.79	184.49
            246 121
            278.76 214];

objPoint = [323 59  		
            300.118 112.8	
            183 109 		
            278.5 55.0
            292 166.4];


objPoint = objPoint;
Hobj2img = calH(objPoint,imgPoint);

pointNum = size(imgPoint,1);
for i = 1:pointNum
    op = [objPoint(i,1) objPoint(i,2) 1]';
    ip = Hobj2img * op;
    ip = ip/ip(3);
    fprintf('%f, %f\n', ip(1),ip(2));
end

width = 448;
height = 200;
correctImg = zeros(height, width, 3);
ori_width = size(img,2);
ori_height = size(img,1);
for iRow = 1:height
    for iCol = 1:width
       Xcorr = [iCol,iRow, 1]';
       X = Hobj2img*Xcorr;
       X = int32(X/X(3));       
       row = X(2);
       col = X(1);
       if(row < 1 | row > ori_height | col< 1 |col > ori_width)
%           correctImg(iRow, iCol,4) = 0; 
           continue;
       end
       r = img(row,col,1);
       g = img(row,col,2);
       b = img(row,col,3);
       
       correctImg(iRow, iCol,1) = r; 
       correctImg(iRow, iCol,2) = g; 
       correctImg(iRow, iCol,3) = b; 
%       correctImg(iRow, iCol,4) = 255; 
    end    
end
correctImg = uint8(correctImg);
imshow(correctImg);
imwrite(correctImg, 'correctImg.jpg');
fprintf('img correct Done!\n');


 
 