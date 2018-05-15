function Skeletonization_square( S,tmp,threshold)

[density_H,density_W]=size(S);

[ridges,ridgemap]=Compute_ridges(S,density_H,density_W,tmp,threshold);


figure;imshow((ridgemap+tmp)/2);title('skeleton_map');
end

