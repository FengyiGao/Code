
% % -------Gao fengyi-----% %
% 这里是为了在热曲面显示成3d脊点
%  objmap 最外圈是包含着 boundarymap
function [ ridgepoint_all,ridgepoint1 ,R_map1 ,ridgepoint2, R_map2] = Correct_direction_test3d( direct_max,direct_min,K,P,fx,fy,fxx,fxy,fyy,fxxx,fyyy,fxxy,fxyy,H,W,objmap,map3d,boundarymap)
% input
%   direct_max direct_min 分别表示最大主方向和最小主方向
%   K 表示n*2的特征值
%   P 表示2*2*n的特征向量
%   H /W表示图像的 高 /宽
%   其余的均表示偏导数
%  boundarymap是边界图 边界处为1 非边界为0
theta=1e-8; %用它来判断两个点的异号
ridgepoint_all=[];
ridgepoint1=[];
ridgepoint2=[];
r_test=[];
R_map1=zeros(H,W);
R_map2=zeros(H,W);
disp('correct direction and compute ridge points');
for i=1:H-1     %先处理横向的边
   for j=1:W-1
       
       ind1=sub2ind([H,W],i,j);
       ind2=sub2ind([H,W],i,j+1);
       if  objmap(i,j)==1 && objmap(i,j+1)==1 && ~(boundarymap(i,j)==1 && boundarymap(i,j+1)==1)
           cor=[i,j,map3d(i,j);i,j+1,map3d(i,j+1)]; % 记录当前处理的两个点的平面坐标 其Z值对应于R

           d1=direct_max(ind1,:);
           d2=direct_max(ind2,:);
           d1_min=direct_min(ind1,:);
           d2_min=direct_min(ind2,:);
           fx_=[fx(ind1);fx(ind2)];
           fy_=[fy(ind1);fy(ind2)];
           fxx_=[fxx(ind1);fxx(ind2)];
           fxy_=[fxy(ind1);fxy(ind2)];
           fyy_=[fyy(ind1);fyy(ind2)];
           fxxx_=[fxxx(ind1);fxxx(ind2)];
           fyyy_=[fyyy(ind1);fyyy(ind2)];
           fxxy_=[fxxy(ind1);fxxy(ind2)];
           fxyy_=[fxyy(ind1);fxyy(ind2)];

           k_two_max=[K(ind1,1);K(ind2,1)]; 
           k_two_min=[K(ind1,2);K(ind2,2)];
           % t=P(:,1,ind1);% test
           if(d1*d2'<0 )
              p_two_max=[P(:,1,ind1),-P(:,1,ind2)]; %最大主曲率对应的
           else
               p_two_max=[P(:,1,ind1),P(:,1,ind2)];
           end

           if (d1_min*d2_min'<0)
               p_two_min=[P(:,2,ind1),-P(:,2,ind2)];
           else
               p_two_min=[P(:,2,ind1),P(:,2,ind2)];
           end
           R_max=compute_R(k_two_max,p_two_max,fx_,fy_,fxx_,fxy_,fyy_,fxxx_,fyyy_,fxxy_,fxyy_);
           R_map1(ind1)=R_max(1);
           R_map1(ind2)=R_max(2);
           [point,r1]=compute_ridgepoint(cor,R_max);
 
           R_min=compute_R(k_two_min,p_two_min,fx_,fy_,fxx_,fxy_,fyy_,fxxx_,fyyy_,fxxy_,fxyy_);
           R_map2(ind1)=R_min(1);
           R_map2(ind2)=R_min(2);
           point_min=compute_ridgepoint(cor,R_min);
           if ~isempty(point)
               ridgepoint_all=[ridgepoint_all;point];
                 if xor(boundarymap(cor(1,1),cor(1,2))==1,boundarymap(cor(2,1),cor(2,2))==1)% xor(boundarymap(cor(1,1),cor(1,2)),boundarymap(cor(2,1),cor(2,2)))
                      if boundarymap(cor(1,1),cor(1,2))==1
                          if sqrt(sum((point(:,1:2)-cor(1,1:2)).^2))<0.5
                              point=[];r1=[];
                          end
                      else
                          if sqrt(sum((point(:,1:2)-cor(2,1:2)).^2))<0.5
                              point=[];r1=[];
                          end
                      end

                 end
              ridgepoint1=[ridgepoint1;point]; 
              r_test=[r_test;r1'];
           end
           if ~isempty(point_min)
              ridgepoint2=[ridgepoint2;point_min] ;
           end
       end
       
       %下面计算纵向的边
       if objmap(i,j)==1 && objmap(i+1,j)==1 && ~(boundarymap(i,j)==1 && boundarymap(i+1,j)==1)
%            if (i==99 && j==114) 
%                disp('');
%            end
           ind1v=sub2ind([H,W],i,j);
           ind2v=sub2ind([H,W],i+1,j);
           corv=[i,j,map3d(i,j);i+1,j,map3d(i+1,j)]; % 记录当前处理的两个点的平面坐标 其Z值对应于R

           d1=direct_max(ind1v,:);
           d2=direct_max(ind2v,:);
           d1_min=direct_min(ind1v,:);
           d2_min=direct_min(ind1v,:);

           fx_=[fx(ind1v);fx(ind2v)];
           fy_=[fy(ind1v);fy(ind2v)];
           fxx_=[fxx(ind1v);fxx(ind2v)];
           fxy_=[fxy(ind1v);fxy(ind2v)];
           fyy_=[fyy(ind1v);fyy(ind2v)];
           fxxx_=[fxxx(ind1v);fxxx(ind2v)];
           fyyy_=[fyyy(ind1v);fyyy(ind2v)];
           fxxy_=[fxxy(ind1v);fxxy(ind2v)];
           fxyy_=[fxyy(ind1v);fxyy(ind2v)];

           k_two_max=[K(ind1v,1);K(ind2v,1)]; 
           k_two_min=[K(ind1v,2);K(ind2v,2)];
           if(d1*d2'<0 )
              p_two_max=[P(:,1,ind1v),-P(:,1,ind2v)]; %最大主曲率对应的
           else
               p_two_max=[P(:,1,ind1v),P(:,1,ind2v)];
           end
           if(d1_min*d2_min'<0)
               p_two_min=[P(:,2,ind1v),-P(:,2,ind2v)];
           else
               p_two_min=[P(:,2,ind1v),P(:,2,ind2v)];
           end
           R=compute_R(k_two_max,p_two_max,fx_,fy_,fxx_,fxy_,fyy_,fxxx_,fyyy_,fxxy_,fxyy_);
           R_map1(ind1v)=R(1);
           R_map1(ind2v)=R(2);
           [point,r2]=compute_ridgepoint(corv,R);
           
           if ~isempty(point)
               ridgepoint_all=[ridgepoint_all;point];
               if xor( boundarymap(corv(1,1),corv(1,2))==1,boundarymap(corv(2,1),corv(2,2))==1)
                  if boundarymap(corv(1,1),corv(1,2))==1
                      if sqrt(sum((point(:,1:2)-corv(1,1:2)).^2))<0.5
                          point=[];r2=[];
                      end
                  else
                      if sqrt(sum((point(:,1:2)-corv(2,1:2)).^2))<0.5
                          point=[];r2=[];
                      end
                  end

               end

                  ridgepoint1=[ridgepoint1;point]; 
                  r_test=[r_test;r2'];
          end

           R_min=compute_R(k_two_min,p_two_min,fx_,fy_,fxx_,fxy_,fyy_,fxxx_,fyyy_,fxxy_,fxyy_);
           R_map2(ind1v)=R_min(1);
           R_map2(ind2v)=R_min(2);
           point_min=compute_ridgepoint(corv,R_min);
           if ~isempty(point_min)
               ridgepoint2=[ridgepoint2;point_min];
           end
  
       
       end
   end 
   
end

%% 测试 代数方法
% figure;plot(ridgepoint1(:,1),ridgepoint1(:,2),'r.');title('all ridges');
% r_new=sort(r_test,2);
% rr1=find(r_new(:,1)<-theta & r_new(:,2)>theta);
% wrong_ridge=ridgepoint1(rr1,:);
% figure;plot(wrong_ridge(:,1),wrong_ridge(:,2),'b.');title('wrong ridges');
% hh=1:size(ridgepoint1,1);
% rr2=setdiff(hh',rr1);
% right_ridge=ridgepoint1(rr2,:);
% figure;plot(right_ridge(:,1),right_ridge(:,2),'g.');title('right ridges');
% figure;
% for i=1:size(ridgepoint1,1)  
%     plot3([ridgepoint1(i,1);ridgepoint1(i,1)],[ridgepoint1(i,2);ridgepoint1(i,2)],[r_new(i,1);r_new(i,2)],'r-');
%     hold on; 
% end
%  hold on;
%  plot3(ridgepoint1(:,1),ridgepoint1(:,2),r_new(:,2),'b.'); % 正R
%  hold on;
%  plot3(ridgepoint1(:,1),ridgepoint1(:,2),r_new(:,1),'g.') % 负R

end

function R=compute_R(K_two,P_two,fx,fy,fxx,fxy,fyy,fxxx,fyyy,fxxy,fxyy)

p1=P_two(1,:);
p2=P_two(2,:);
p1=p1(:);
p2=p2(:);
t1=p1.^3;
t2=(p1.^2).*p2;
t3=p1.*(p2.^2);
t4=p2.^3;
t5=p1.^2;
t6=p1.*p2;
t7=p2.^2;
t8=p1;
t9=p2;

k=K_two;

R=(t1.*fxxx+(3*t2).*fxxy+(3*t3).*fxyy+t4.*fyyy)-3*sqrt(ones(size(fx))...
   + fx.*fx + fy.*fy).* (t5.*fxx+(2*t6).*fxy + t7.*fyy) .*( t8.*fx+t9 .*fy ).*k ;

end

function [point,rtmp]=compute_ridgepoint(cor,R) 
% 利用辛老师说的分比来求0点
point=[];
rtmp=[];
if (R(1)*R(2)<0 )% R(1)*R(2)<=0 %(R(1)*R(2)<0 || R(1)*R(2)==0)
    rate=abs(R(1))./abs(R(1)-R(2));
    a=[cor(1,:)];
    b=[cor(2,:)];
    ab=b-a;
    point=rate*ab+a;
    rtmp=R;
end

end




