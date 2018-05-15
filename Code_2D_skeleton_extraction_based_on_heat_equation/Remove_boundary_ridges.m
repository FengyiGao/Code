function [ new_ridges ] = Remove_boundary_ridges( all_ridges,boundary,threshold )
%REMOVE_BOUNDARY_RIDGES 此处显示有关此函数的摘要
%   此处显示详细说明

[H,h]=size(all_ridges);
if h>2
    new_ridges=all_ridges(:,1:2);
else
    new_ridges=all_ridges;
end
s=true(H,1);
for i=1:size(new_ridges,1)
   cur_ridge=new_ridges(i,:) ;
   d=sqrt(sum((repmat(cur_ridge,size(boundary,1),1)-boundary).^2,2));
   if min(d)<threshold * 0.99 %这里设置了一个容错概率 %|| min(d)==threshold
      s(i)=false;
   end
end
new_ridges=new_ridges(s,:);

end

