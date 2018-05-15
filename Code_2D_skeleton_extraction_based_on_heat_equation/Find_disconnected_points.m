function [ new_ridgemap ] = Find_disconnected_points( new_ridges,H,W,threshold ,map3d,objmap)
%FIND_DISCONNECTED_POINTS 此处显示有关此函数的摘要
%   threshold 用来截断一些小的分支对应的脊点
% sign=true(size(new_ridges,1),1);
ridgemap=zeros(H,W);

%% 用四舍五入 这样得到的基本是单像素的骨架图
inter_ridge=round(new_ridges);

if isempty(inter_ridge)
    disp('none ridges');
else
    inter_ridge_ind=sub2ind([H,W],inter_ridge(:,1),inter_ridge(:,2));
try
   ridgemap(inter_ridge_ind)=1;
catch
   disp('wrong in find disconnected points') ;
end
%%  
        disconnected_uv=[];
        disconnected=[];
        neigh8=[0,-1;0,1;1,0;1,1;1,-1;-1,0;-1,1;-1,-1];

        [L,num]=bwlabel(ridgemap,8);
           for j=1:num
               if num==1
                   break;
               end
              cur_j=find(L==j);
              if size(cur_j,1)<=threshold
                 ridgemap(cur_j)=0; 
              end
           end
        %%
        thin_ridgemap=bwmorph(ridgemap,'thin');
        thin_inter_ridge_ind=find(thin_ridgemap);
        [thin_u,thin_v]=ind2sub(size(thin_ridgemap),thin_inter_ridge_ind);
        for i=1:size(thin_inter_ridge_ind,1)
            cur=[thin_u(i),thin_v(i)];
            cur_neigh=repmat(cur,size(neigh8,1),1)+neigh8;
            s=cur_neigh(:,1)<1 |cur_neigh(:,1)>H |cur_neigh(:,2)<1 |cur_neigh(:,2)>W ;
            cur_neigh(s,:)=[];
            cur_neigh_ind=sub2ind([H,W],cur_neigh(:,1),cur_neigh(:,2));
            iii=find(thin_ridgemap(cur_neigh_ind));
            if size(iii,1)<=1
                disconnected_uv=[disconnected_uv;cur];
            end
                
        end
         old_ind=find(ridgemap);
         

%%
  [new_ridgemap]= trace(disconnected,disconnected_uv,map3d,ridgemap);

end
end

function [new_ridgemap]=trace(disconnected,disconnected_uv,map3d,ridgemap)

    % 从断头开始trace
    pixel=true;
    new_ridgemap=ridgemap;
    
    [H,W]=size(ridgemap);
    neigh8=[-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1];
    ridgemap_l=logical(ridgemap);
    if pixel % 使用像素trace 8邻域比较大小
        for i=1:size(disconnected_uv,1)
            cur=disconnected_uv(i,:);
            while(1)
            cur_neigh=repmat(cur,size(neigh8,1),1)+neigh8;
            ind=cur_neigh(:,1)<1 | cur_neigh(:,1)>H |cur_neigh(:,2)<1 |cur_neigh(:,2)>W;
            if any(ind)
                cur_neigh(ind,:)=[];
            end
            cur_neigh_ind=sub2ind(size(ridgemap),cur_neigh(:,1),cur_neigh(:,2));
            cur_neigh_value=map3d(cur_neigh_ind);
            [~,max_ind]=max(cur_neigh_value);
             if ridgemap_l(cur_neigh_ind(max_ind))
                break;
            end  %四周有标记点 则为相撞
            new_ridgemap(cur_neigh_ind(max_ind))=1;
            ridgemap_l(cur_neigh_ind(max_ind))=true;
            cur=cur_neigh(max_ind,:);
            end
        end
    
        
    end

end

