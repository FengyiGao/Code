% 重新整理一下计算脊点的部分 
function [ new_ridges,new_ridgemap] = Compute_ridges( map3d,H,W,objmap,threshold_ridge)

% 利用公式计算骨架点位置
% 输入V 的点的xy坐标都是整数值 z是对应于热方程计算出来的值
% map表示经过膨胀后的新的01图像 也可能未膨胀
% map_b表示的是分辨率增加后边界图像 边界上的值为1 其他值为0
% mainarea是经过填充后的显著性区域的单标号 和V[:,1:2]应该是一一对应 只是V多了外圈

%% 一些初始化
bounda=bwboundaries(objmap);
boundary=[];
for ab=1:size(bounda,1)
    boundary=[boundary;bounda{ab,1}];
end
b_ind=sub2ind([H,W],boundary(:,1),boundary(:,2));
map_b=zeros(H,W);
map_b(b_ind)=1;

%%
sub_x=boundary(:,1);
sub_y=boundary(:,2);
% %尝试把x y反过来试试 前后列相减得到fx 行相减fy

%% 刚注释掉以下的向前向后差分 这里是中心差分
% 这边要注意看x,y的方向 fx-表示图像坐标系中的列之差 fy-表示图像坐标系中的行之差
[ fx,fy,fxx,fyy,fxy,fyxx,fxyy,fxxx,fyyy ] =Central_diff(map3d);
fx_=fx(:);
fy_=fy(:);
fxx_=fxx(:);
fyy_=fyy(:);
fxy_=fxy(:);
fxxy_=fyxx(:);
fxyy_=fxyy(:);
fxxx_=fxxx(:);
fyyy_=fyyy(:);

iidx=1;
% iidx=sub2ind(size(map3d),199,142);
%% 开始计算脊点的一些参数 p1 p2 kp
% P 是2*2*n  K 是n*2 大的特征值在第一列 大的特征值对应的特征向量在第一列

[P,K,direct_max,direct_min,normal,other]=compute_ridge_parameters(fx_,fy_,fxx_,fxy_,fyy_,iidx);
%% 矫正相邻两点的主方向 若>90°，p1 p2取反 ，同时直接计算两点间是否存在脊点
[all_ridge,ridge_points_max,R_max,ridge_points_min,R_min]=Correct_direction_test3d(direct_max,direct_min,K,P,fx_,fy_,fxx_,fxy_,fyy_,fxxx_,fyyy_,fxxy_,fxyy_,H,W,objmap,map3d,map_b);

new_ridges=Remove_boundary_ridges(ridge_points_max,boundary,threshold_ridge);

figure;
plot(new_ridges(:,1),new_ridges(:,2),'r.');
hold on
plot(sub_x,sub_y,'b.');
axis equal;axis off
title('new threshold');

[new_ridgemap]=Find_disconnected_points(new_ridges(:,1:2),H,W,4,-map3d,objmap);

end




