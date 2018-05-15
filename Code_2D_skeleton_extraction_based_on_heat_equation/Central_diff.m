function [ fx,fy,fxx,fyy,fxy,fyxx,fxyy,fxxx,fyyy ] = Central_diff( map )
%CENTRAL_DIFF 利用中心差分求一阶导数
% 后面在处理的时候是直接在 xoy坐标中进行 所以现在 1<=x<row 1<=y<=col 假设uov坐标原点和xoy原点重合
step=1;
map_xoy=map; %(i,j) 像素坐标，不是笛卡尔坐标
iandjcut1=[map_xoy(:,1),map_xoy(:,1:end-1)];% (i,j-1)
iandjadd1=[map_xoy(:,2:end),map_xoy(:,end)];% (i,j+1)

iadd1andj=[map_xoy(2:end,:);map_xoy(end,:)]; %(i+1,j)
icut1andj=[map_xoy(1,:);map_xoy(1:end-1,:)]; %(i-1,j)

icut1andjadd1=[map_xoy(1,:);map_xoy(1:end-1,2:end),map_xoy(1:end-1,end)]; %(i-1,j+1)
iadd1andjadd1=[map_xoy(2:end,2:end),map_xoy(2:end,end);map_xoy(end,:)]; % (i+1,j+1)

icut1andjcut1=[map_xoy(1,:);map_xoy(1:end-1,1),map_xoy(1:end-1,1:end-1)]; %(i-1,j-1)
iadd1andjcut1=[map_xoy(2:end,1),map_xoy(2:end,1:end-1);map_xoy(end,:)]; %(i+1,j-1)

iadd2andj=[map_xoy(3:end,:);map_xoy(end-1:end,:)]; %(i+2,j)
icut2andj=[map_xoy(1:2,:);map_xoy(1:end-2,:)];%(i-2,j)
iandjcut2=[map_xoy(:,1:2),map_xoy(:,1:end-2)];%(i,j-2)
iandjadd2=[map_xoy(:,3:end),map_xoy(:,end-1:end)];%(i,j+2)

fy=(iadd1andj-icut1andj)/2*step;  %fy(i,j)=(f(i+1,j)-f(i-1,j))/2
fx=(iandjadd1-iandjcut1)/2*step;  %fx(i,j)=(f(i,j+1)-f(i,j-1))/2                                                                    

fyy=iadd1andj-2*map_xoy+icut1andj; %fyy(i,j)=f(i+1,j)+f(i-1,j)-2*f(i,j)
fxx=iandjadd1-2*map_xoy+iandjcut1; %fxx(i,j)=f(i,j+1)+f(i,j-1)-2*f(i,j)

fxy=(iadd1andjadd1-icut1andjadd1-iadd1andjcut1+icut1andjcut1)/4; %fxy(i,j)=(f(i+1,j+1)+f(i-1,j-1)-f(i-1,j+1)-f(i+1,j-1))/4
fxxx=(iandjadd2-2*iandjadd1+2*iandjcut1-iandjcut2)/2; %fxxx(i,j)=f(i,j+2)/2-f(i,j-2)/2+f(i,j-1)-f(i,j+1)
fyyy=(iadd2andj-2*iadd1andj+2*icut1andj-icut2andj)/2; %fyyy(i,j)=f(i-1,j)-f(i+1,j)-f(i-2,j)/2+f(i+2,j)/2
fyxx=(iadd1andjadd1-2*iadd1andj+iadd1andjcut1-icut1andjadd1+2*icut1andj-icut1andjcut1)/2; % fyxx(i,j)=-f(i,j+1)+f(i-1,j+1)/2+f(I+1,j+1)/2-f(i-1,j-1)/2+f(i,j-1)-f(i+1,j-1)/2
fxyy=(iadd1andjadd1-2*iandjadd1+icut1andjadd1-iadd1andjcut1+2*iandjcut1-icut1andjcut1)/2; 
%假设 横着是dx 纵着是dy



end

