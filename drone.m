%% Init
% clc;clear;close;
imgNm = 'map_in.png';
img = imread(imgNm);

% figure(1);
% imshow(img);

hsv = rgb2hsv(img);
[x,y,~] = size(img);

%% Retrieve the region
imgin = zeros(size(img));
for i=1:x
    for j=1:y
        t = hsv(i,j,:);
        if t(2) < 0.01 && t(3) > 0.9
            imgin(i,j,1) = 1;
        end
    end
end

% figure(2);
% imshow(imgin);

%% Load the road map
rdmpNm = 'rdmp.png';
rdmp = imread(rdmpNm);

% figure(3);
% imshow(rdmp);

% Fix the blure and normalize
[x,y,~] = size(rdmp);
hsv = rgb2hsv(rdmp);
for i=1:x
    for j=1:y
        t = hsv(i,j,:);
        if t(2) > 0.5 && t(3) > 0.5
            rdmp(i,j,1) = 255;
            rdmp(i,j,2) = 0;
            rdmp(i,j,3) = 0;
        else
            rdmp(i,j,1) = 0;
            rdmp(i,j,2) = 0;
            rdmp(i,j,3) = 0;
        end
        
%         if rdmp(i,j,2) > 50 && rdmp(i,j,2) > 50
%         if rdmp(i,j,1) < 40
%             rdmp(i,j,1) = 0;
%         else
%             rdmp(i,j,1) = 255;
%         end
    end
end

% figure(4);
% imshow(rdmp);
rdmp = im2single(rdmp);

%% Calculate the value of several points
vm = 0;
xi = 0;yi = 0;
for i=1:x
    if mod(i,50)==0
        disp(['calculating...',num2str(i),'/679']);
    end
    for j=1:y
        if imgin(i,j,1) == 1
            vans = calQ(i,j,rdmp);
            if vm <= vans
                vm = vans;
                xi = i;yi = j;
            end
        end
    end
end

%% Visualization
% vm
% xi
% yi
wid = 2;
lth = 17;
% imgout = imread('map.jpg');
% imgout = imread('map.png');
imgout = imread('combine.png');
tr = 255;
tg = 0;
tb = 170;
for i = -lth:1:lth
    imgout(xi+i,yi-wid:yi+wid,1) = tr;
    imgout(xi+i,yi-wid:yi+wid,2) = tg;
    imgout(xi+i,yi-wid:yi+wid,3) = tb;
end
for i = -lth:1:lth
    imgout(xi-wid:xi+wid,yi+i,1) = tr;
    imgout(xi-wid:xi+wid,yi+i,2) = tg;
    imgout(xi-wid:xi+wid,yi+i,3) = tb;
end
figure(16);
imshow(imgout);
imwrite(imgout,'imgout.png');

%% Utility function
function vans = calQ(cx,cy,rdmp)
    R = 100;
    R2 = R^2;
    vans = 0;
    il = max(1,cx-R);
    ir = min(679, cx+R);
    jl = max(1,cy-R);
    jr = min(1531,cy+R);
    for i=il:ir
        for j=jl:jr
            if (i-cx)^2 + (j-cy)^2 < R2
                vans = vans + double(rdmp(i,j,1));
            end
        end
    end
end

