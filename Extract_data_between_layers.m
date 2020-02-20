function [structure_withline,structure_Flattening,floe_Flattening,retinalLayers_mat]=Extract_data_between_layers(structure,floe,num1,num2,retinalLayers_mat)
% 根据oct结构图的分层结果把血流图的num1层和num2层之间的图像保存下来。
%并且根据第四层拉平
%输入：
%      structure：结构图矩阵
%      floe：血流图矩阵
%      num1：上层序号 num1==1时，num1层以上不置0
%      num2：下层序号 num2==8时，num8层以下不置0
%输出：
%      floe_extract：保存的分层
% Author : Ming, 9/12/2019
    structure=mat2gray(structure);
    structure=structure./max(structure(:));

    if nargin < 5 % 如果输入参数中没有retinalLayers_mat，那么重新计算分层
         retinalLayers_mat=Layers_location(imgaussfilt3(structure,3));%三维中值滤波
    end
    [height,width,depth]=size(structure); 
    structure_withline(:,:,:,1)=structure.*1;
    structure_withline(:,:,:,2)=structure.*1;
    structure_withline(:,:,:,3)=structure.*1;
    structure_Flattening=structure;%拉平的结构图
    floe_Flattening=floe;%拉平的血流图
    [height_retinal,width_retinal,depth_retinal]=size(retinalLayers_mat); 
    colorarr=colormap('jet'); 
    colorarr=colorarr(64:-8:1,:);%使用到的colorbar
    colorarr1=colorarr; 
    %1 2 3 4     
    %5 6 7 8
    colorarr1(1,:)=[0,0,0];
    colorarr1(2,:)=[1,0,0];
    colorarr1(3,:)=[0,1,0];
    colorarr1(4,:)=[0,0.9,0];
    colorarr1(5,:)=[0,0.8,0];
    colorarr1(6,:)=[0,0.7,0];
    colorarr1(7,:)=[1,0,0];
    colorarr1(8,:)=[1,0,0];
    colorarr=colorarr1;
    %逐frame 根据分层信息进行处理
    for kk=1:depth
        drawline=structure_withline(:,:,kk,:);
        structure_image=structure_Flattening(:,:,kk);%结构图
        floe_image=floe_Flattening(:,:,kk);%血流图
        
        for jj = 1 : width
            %将每一个A-line的分层保存在location
            for iii=1:height_retinal
                location(iii)=retinalLayers_mat(iii,jj,kk);
            end
            %根据location中保存的分层信息，对drawline赋予不同颜色
            for iii=1:height_retinal-1
                colora = colorarr(iii,:);
                %为分层不同的颜色
%                 drawline(location(iii):location(iii+1), jj,1)=structure_image(location(iii):location(iii+1), jj).*colora(1);%[1,0,0] 
%                 drawline(location(iii):location(iii+1), jj,2)=structure_image(location(iii):location(iii+1), jj).*colora(2);%[1,0,0] 
%                 drawline(location(iii):location(iii+1), jj,3)=structure_image(location(iii):location(iii+1), jj).*colora(3);%[1,0,0] 
                %分层分界线不同颜色
                drawline(location(iii):location(iii)+1, jj,1)=colora(1);%[1,0,0] 
                drawline(location(iii):location(iii)+1, jj,2)=colora(2);%[1,0,0] 
                drawline(location(iii):location(iii)+1, jj,3)=colora(3);%[1,0,0] 
            end
            %血流
            im=drawline(location(3):location(7), jj,3);             % channel 3
            im(floe_image(location(3):location(7),jj)>0.08)=1;
            drawline(location(3):location(7), jj,1)=im;

            if (num1~=1)
                %num1层以上置零
                location=retinalLayers_mat(num1,jj,kk);
%                 structure_image(1:location+7, jj)=0;
                if num1==4
                    floe_image(1:location-0, jj)=0;
                else
                    floe_image(1:location, jj)=0;
                end
            end
            if (num2~=height_retinal)
                %num2层以下置零
                location=retinalLayers_mat(num2,jj,kk);
%                 structure_image(location-10:height, jj)=0;
                floe_image(location-1:height, jj)=0;
            end
            %根据第四层拉平
            location=retinalLayers_mat(height_retinal-2,jj,kk);
            drawline(:, jj,1) = circshift(drawline(:, jj,1),  round(height/2)- location);
            drawline(:, jj,2) = circshift(drawline(:, jj,2),  round(height/2)- location);
            drawline(:, jj,3) = circshift(drawline(:, jj,3),  round(height/2)- location);

            structure_image(:, jj) = circshift(structure_image(:, jj),  round(height/2)- location);
            floe_image(:, jj) = circshift(floe_image(:, jj),  round(height/2)- location);
        end
        structure_withline(:,:,kk,:)=drawline;
        structure_Flattening(:,:,kk)=structure_image;
        floe_Flattening(:,:,kk)=floe_image;
        imshow(floe_image)
    end
    
    
end