% Layere_Division
%对有中央大凹陷的眼底图像分层，程序有待改善

% clear
% structure_image=imread("/Users/ming/Documents/data/DataPicture/HD5-Line2/03355_Meditco_OCT/03355_Meditco_OCT_003.jpg");
% structure_image=imread("/Users/ming/Documents/data/DataPicture/HD5-Line1/01152_Meditco_OCT/01152_Meditco_OCT_003.jpg");
% structure_image=imread("/Users/ming/Documents/data/DataPicture/HD5-Line2/03359_Meditco_OCT/03359_Meditco_OCT_001.jpg");
% structure_image=rgb2gray(structure_image);

[height,width]=size(structure_image);
img = structure_image;
structure_img = img;

img=double(img)./255.0;                                                                                
img(img<1.0*mean(mean(img)))=0;
% img=imbothat(img,strel('sphere',50));
% img=img-imtophat(img,strel('sphere',5));
% imshow(img)
% img=sqrt(img);
% img=img.*img;
% img=medfilt2(img,[3,3]);%中值滤波
img = imfilter(img, fspecial('gaussian', [15 5], 3), 'symmetric');
rawimg =img;

%% 找到凹陷
[~, params] = getRetinalLayers(rawimg);
img = imfilter(img,fspecial('gaussian',params.filter0Params(1:2),params.filter0Params(3)),'replicate');        
[params.adjMatrixW, params.adjMatrixMW, params.adjMA, params.adjMB, params.adjMW, params.adjMmW, imgNew] = getAdjacencyMatrix(img);
retinalLayerSegmentationOrder = {'roughILMandISOS' 'ilm'};
retinalLayers = [];
distance=[0,1,0,0,0];%1：记录凹陷处偏移两端点连线的距离 2：凹陷处坐标 3：边缘每一个点到两端点连线的平均距离
for layerInd = 1:numel(retinalLayerSegmentationOrder)
    %找到了每一层分界线
    [retinalLayers, ~,distance] = getRetinalLayersCore(retinalLayerSegmentationOrder{layerInd},imgNew,params,retinalLayers,0,distance);
end

if distance(1)>10 && abs(distance(4)-distance(5))>200
    img_1=rawimg;
    img_1(:,distance(4):distance(5))=[];
    imshow(img_1);
    [retinalLayers1, params1] = getRetinalLayers(img_1);
    img_2=rawimg(:,distance(4):distance(5));
    [retinalLayers2, params2] = getRetinalLayers(img_2);
    img_3=rawimg(:,distance(5):width);
    [retinalLayers3, params3] = getRetinalLayers(img_3);
    
    layersToPlot = {'ilm' 'isos' 'rpe' 'inlopl' 'nflgcl' 'iplinl' 'oplonl'};% 'rpeSmooth'}; %
    retinalLayers=retinalLayers1;
    for k = 1:numel(layersToPlot)
        matchedLayers = strcmpi(layersToPlot{k},{retinalLayers(:).name});
        layerToPlotInd = find(matchedLayers == 1);
        index1=length(retinalLayers1(layerToPlotInd).pathY());
        index2=length(retinalLayers2(layerToPlotInd).pathY());
        index3=length(retinalLayers3(layerToPlotInd).pathY());
        
        retinalLayers(layerToPlotInd).pathY(1:index1)=retinalLayers1(layerToPlotInd).pathY();
%         retinalLayers(layerToPlotInd).pathY(1:index1)=retinalLayers1(layerToPlotInd).pathX();
        retinalLayers(layerToPlotInd).pathY(index1+1:index1+index2)=retinalLayers2(layerToPlotInd).pathY();
%         retinalLayers(layerToPlotInd).pathY(index1+1:index1+index2)=retinalLayers2(layerToPlotInd).pathX();
        retinalLayers(layerToPlotInd).pathY(index1+index2+1:index1+index2+index3)=retinalLayers3(layerToPlotInd).pathY();
%         retinalLayers(layerToPlotInd).pathY(index1+index2+1:index1+index2+index3)=retinalLayers3(layerToPlotInd).pathY();
    end
else
    [retinalLayers, params] = getRetinalLayers(rawimg);
    % title("Image Segmentation");
end


% plot oct image and the obtained retinal layers.
isPlot = 1;
if isPlot,

    imagesc(img);
    axis image; colormap('gray'); hold on; drawnow;
    layersToPlot = {'ilm' 'isos' 'rpe' 'inlopl' 'nflgcl' 'iplinl' 'oplonl'};% 'rpeSmooth'}; %

%     layersToPlot = {'ilm' 'isos' 'rpe' 'nflgcl' 'iplinl'};% 'rpeSmooth'}; %
    hOffset =       [40    0      40    0        0        40       -40      -40]; % for displaying text
    for k = 1:numel(layersToPlot)

        matchedLayers = strcmpi(layersToPlot{k},{retinalLayers(:).name});
        layerToPlotInd = find(matchedLayers == 1);

        if ~isempty(retinalLayers(layerToPlotInd).pathX)         
            colora = params1.colorarr(k,:);
%             retinalLayers(layerToPlotInd).pathX = smoothdata(retinalLayers(layerToPlotInd).pathX,'gaussian',25);
            retinalLayers(layerToPlotInd).pathY=retinalLayers(layerToPlotInd).pathY+100;
            
            plot(retinalLayers(layerToPlotInd).pathY,retinalLayers(layerToPlotInd).pathX-1,'-','color',colora,'linewidth',1.5);
            plotInd = round(numel(retinalLayers(layerToPlotInd).pathX)/2);            
            text(retinalLayers(layerToPlotInd).pathY(plotInd)+hOffset(k),retinalLayers(layerToPlotInd).pathX(plotInd)+params1.txtOffset,retinalLayers(layerToPlotInd).name,'color',colora,'linewidth',2);            
            drawnow;
        end % of if ~isempty            

    end % of k
    hold off;

end % of isPlot        
