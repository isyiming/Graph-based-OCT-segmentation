function retinalLayers_mat2=Layers_location(structure)
% ��oct�ṹͼ�ķֲ�������������
%���룺
%      structure���ṹͼ����
%�����
%      retinalLayers_mat������ķֲ�����
% ��� structureΪ512*250*200 �ľ���cro-section��СΪ512*250
% ��ôretinalLayers_mat��СΪ5*250*200����ȡ�Ľṹͼ��һ���������ΪretinalLayers_mat(1,:,:)
% Author : Ming, 9/12/2019

    [height,width,depth]=size(structure); 
    %retinalLayers ���ص���7�� retinalLayers_mat��Ϊ9�㣬��һ�����Ϊ0 �ھŲ����Ϊheight
    retinalLayers_mat=zeros(9,width,depth);
    retinalLayers_mat(1,:,:)=1;
    retinalLayers_mat(9,:,:)=height;
    [height_retinal,width_retinal,depth_retinal]=size(structure); 
    
    %��frame ���ݷֲ���Ϣ���д���
    for kk=1:depth
        img=structure(:,:,kk);%����߽�Ľṹͼ
        [retinalLayers, params1] = getRetinalLayers(img);%�ֲ�
        %���õ��ֲ���Ϣ
        sz_retinal=size(retinalLayers);
        for ii=1:sz_retinal(2)
            pathX=retinalLayers(ii).pathX;
            pathY=retinalLayers(ii).pathY; 
             %ÿһ��A-line�ķֲ�
            for jj = 1 : width
                location=find(pathY==jj);
                retinalLayers_mat(ii+1,jj,kk)=pathX(location(1));%��Ϊ��һ�������Ϊ1��
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
%         retinalLayers_mat(ii,:,:) = medfilt3(retinalLayers_mat(ii,:,:),[1,3,3]);%��ֵ�˲�
        retinalLayers_mat2(ii,:,:) = imgaussfilt3(retinalLayers_mat2(ii,:,:), 1);%��˹�˲�
        map=reshape(retinalLayers_mat2(ii,:,:),[width,depth]);
        mesh(map)
        hold on;
    end
    retinalLayers_mat2=round(retinalLayers_mat2);
    retinalLayers_mat2(retinalLayers_mat2==0)=1;
    hold off;
end