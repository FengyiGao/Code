% ��������һ�¼��㼹��Ĳ��� 
function [ new_ridges,new_ridgemap] = Compute_ridges( map3d,H,W,objmap,threshold_ridge)

% ���ù�ʽ����Ǽܵ�λ��
% ����V �ĵ��xy���궼������ֵ z�Ƕ�Ӧ���ȷ��̼��������ֵ
% map��ʾ�������ͺ���µ�01ͼ�� Ҳ����δ����
% map_b��ʾ���Ƿֱ������Ӻ�߽�ͼ�� �߽��ϵ�ֵΪ1 ����ֵΪ0
% mainarea�Ǿ������������������ĵ���� ��V[:,1:2]Ӧ����һһ��Ӧ ֻ��V������Ȧ

%% һЩ��ʼ��
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
% %���԰�x y���������� ǰ��������õ�fx �����fy

%% ��ע�͵����µ���ǰ����� ���������Ĳ��
% ���Ҫע�⿴x,y�ķ��� fx-��ʾͼ������ϵ�е���֮�� fy-��ʾͼ������ϵ�е���֮��
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
%% ��ʼ���㼹���һЩ���� p1 p2 kp
% P ��2*2*n  K ��n*2 �������ֵ�ڵ�һ�� �������ֵ��Ӧ�����������ڵ�һ��

[P,K,direct_max,direct_min,normal,other]=compute_ridge_parameters(fx_,fy_,fxx_,fxy_,fyy_,iidx);
%% ������������������� ��>90�㣬p1 p2ȡ�� ��ͬʱֱ�Ӽ���������Ƿ���ڼ���
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




