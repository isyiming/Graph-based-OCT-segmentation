function retinalLayers_mat2=Layers_location(structure)
% 把oct结构图的分层结果保存下来。
%输入：
%      structure：结构图矩阵
%输出：
%      retinalLayers_mat：保存的分层坐标
% 如果 structure为512*250*200 的矩阵，cro-section大小为512*250
% 那么retinalLayers_mat大小为5*250*200，提取的结构图第一层深度坐标为retinalLayers_mat(1,:,:)
% Author : Ming, 9/12/2019

    [height,width,depth]=size(structure); 
    %retinalLayers 返回的是7层 retinalLayers_mat设为9层，第一层深度为0 第九层深度为height
    retinalLayers_mat=zeros(9,width,depth);
    retinalLayers_mat(1,:,:)=1;
    retinalLayers_mat(9,:,:)=height;
    [height_retinal,width_retinal,depth_retinal]=size(structure); 
    
    %逐frame 根据分层信息进行处理
    for kk=1:depth
        img=structure(:,:,kk);%计算边界的结构图
        [retinalLayers, params1] = getRetinalLayers(img);%分层
        %逐层得到分层信息
        sz_retinal=size(retinalLayers);
        for ii=1:sz_retinal(2)
            pathX=retinalLayers(ii).pathX;
            pathY=retinalLayers(ii).pathY; 
             %每一个A-line的分层
            for jj = 1 : width
                location=find(pathY==jj);
                retinalLayers_mat(ii+1,jj,kk)=pathX(location(1));%因为第一层深度设为1了
            end
%             clf;
%             plot(pathY,pathX);
%             hold on;
%             plot(retinalLayers_mat(ii+1,:,kk))
        end
    end
    retinalLayers_mat2=retinalLayers_mat;
    figure(10);
    for ii=1:8
        mat=retinalLayers_mat(ii,:,:);
        mat=reshape(mat,[width,depth]);
%         retinalLayers_mat2(ii,:,:)=MedianFilterWithOriginalImage(mat);
%         retinalLayers_mat(ii,:,:) = medfilt3(retinalLayers_mat(ii,:,:),[1,3,3]);%中值滤波
        retinalLayers_mat2(ii,:,:) = imgaussfilt3(retinalLayers_mat2(ii,:,:), 1);%高斯滤波
        map=reshape(retinalLayers_mat2(ii,:,:),[width,depth]);
        mesh(map)
        hold on;
    end
    retinalLayers_mat2=round(retinalLayers_mat2);
    retinalLayers_mat2(retinalLayers_mat2==0)=1;
    hold off;
end