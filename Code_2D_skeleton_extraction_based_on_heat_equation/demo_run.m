clc
clear
close all
addpath('data');

thresh=2.0616;
name='obj3_2';

pic=imread(strcat(name ,'.png'));
pic=im2bw(pic);
if pic(1,1)
    pic=~pic;
end

S=load(strcat(name,'.mat'));
S=S.s1;


Skeletonization_square(S,pic,thresh);