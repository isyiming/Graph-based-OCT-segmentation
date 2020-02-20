function [retinalLayers, params] = getRetinalLayers(img,params)
if nargin < 1
    display('requires 1 input');
    return;
end

%initialize constants
if nargin < 2        
    
    % resize the image if 1st value set to 'true',
    % with the second value to be the scale.
    params.isResize = [false 0.5];
    
    % parameter for smothing the images.
    params.filter0Params = [15 5 1];
    params.filterParams = [20 20 2];           
        
    % constants used for defining the region for segmentation of individual layer
    params.roughILMandISOS.shrinkScale = 0.2;
    params.roughILMandISOS.offsets = -40:40;    
    params.ilm_0 = 4;
    params.ilm_1 = 4;
    params.isos_0 = 4;
    params.isos_1 = 4;
    params.rpe_0 = 0.05;
    params.rpe_1 = 0.05;
    params.inlopl_0 = 0.1; %   0.4;%
    params.inlopl_1 = 0.3; %   0.5;%  
    params.nflgcl_0 = 0.05;%  0.01;
    params.nflgcl_1 = 0.3; %   0.1;
    params.iplinl_0 = 0.6;
    params.iplinl_1 = 0.2;
    params.oplonl_0 = 0.05;%4;
    params.oplonl_1 = 0.5;%4;    
        
    % parameters for ploting
    params.txtOffset = -7;
    colorarr=colormap('jet'); 
    params.colorarr=colorarr(64:-8:1,:);
    
    % a constant (not used in this function, used in 'octSegmentationGUI.m'.)
    params.smallIncre = 2;    
    
end

%clear up matlab's mind
clear retinalLayers

%get image size
szImg = size(img);

%resize image.
if params.isResize(1)
    img = imresize(img,params.isResize(2),'bilinear');
end

%smooth image with specified kernels
%for denosing
img = imfilter(img,fspecial('gaussian',params.filter0Params(1:2),params.filter0Params(3)),'replicate');        

%for a very smooth image, a "broad stroke" of the image
imgSmo = imfilter(img,fspecial('gaussian',params.filterParams(1:2),params.filterParams(3)),'replicate');

% create adjacency matrices and its elements base on the image.
[params.adjMatrixW, params.adjMatrixMW, params.adjMA, params.adjMB, params.adjMW, params.adjMmW, imgNew] = getAdjacencyMatrix(img);

% % [this is not used as the moment] Create adjacency matrices and its elements based on the smoothed image.
% [params.adjMatrixWSmo, params.adjMatrixMWSmo, params.adjMA, params.adjMB, params.adjMWSmo, params.adjMmWSmo, ~] = getAdjacencyMatrix(imgSmo);

% obtain rough segmentation of the ilm and isos, then find the retinal
% layers in the order of 'retinalLayerSegmentationOrder'
%%vvvvvvvvvvvvvvvDO  NOT  CHANGE BELOW LINE (ORDER OF LAYERS SHALL NOT BE CHANGED)vvvvvvvvvvvvvv%%
% retinalLayerSegmentationOrder = {'roughILMandISOS' 'ilm' 'isos' 'rpe' 'isos_refine' 'inlopl' 'nflgcl' 'inlopl_refine' 'iplinl' 'oplonl'};

retinalLayerSegmentationOrder = {'roughILMandISOS' 'ilm' 'isos' 'rpe' 'isos_refine' 'inlopl' 'nflgcl' 'iplinl' 'iplinl_refine' 'oplonl' 'inlopl_refine' 'oplonl_refine'};
% retinalLayerSegmentationOrder = {'roughILMandISOS' 'ilm' 'isos' 'rpe' 'isos_refine' 'inlopl' 'nflgcl'};
% retinalLayerSegmentationOrder = {'roughILMandISOS' 'ilm' 'isos' };

%%^^^^^^^^^^^^^^^DO  NOT  CHANGE ABOVE LINE (ORDER OF LAYERS SHOULD NOT BE CHANGED)^^^^^^^^^^^^^%%
% segment retinal layers
retinalLayers = [];
distance=[0,1,0];
for layerInd = 1:numel(retinalLayerSegmentationOrder)
    %找到了每一层分界线
        %%ming
    [retinalLayers, ~,distance] = getRetinalLayersCore(retinalLayerSegmentationOrder{layerInd},imgNew,params,retinalLayers,0,distance);
    
    imagesc(img);
    axis image; colormap('gray'); hold on; drawnow;
    for ii=1:length(retinalLayers)
        plot(retinalLayers(ii).pathY,retinalLayers(ii).pathX-1,'-','color','r','linewidth',1.5);
    end
%     x=retinalLayers(1).pathX;
%     y=retinalLayers(1).pathY;
% 
%     size_img=size(imgNew);
%     for ii=1:size_img(2)
%         imgNew()
%     end
    
end

%delete elements of the adjacency matrices prior function exit to save memory
toBeDeleted = {'adjMatrixWSmo' 'adjMatrixMWSmo' 'adjMWSmo' 'adjMmWSmo'  'adjMW' 'adjMmW' 'adjMatrixW' 'adjMatrixMW' 'adjMA' 'adjMB'};
for delInd = 1:numel(toBeDeleted)
    params.(toBeDeleted{delInd}) = [];
end


for kk=1:length(retinalLayers)
    for n = 1 : szImg(2)
        location=find(retinalLayers(kk).pathY==n);
        if length(location)>1
            for ii=2:length(location)
                retinalLayers(kk).pathY(location(ii))=0;
                retinalLayers(kk).pathX(location(ii))=0;
            end
        end
    end
end
retinalLayers1=retinalLayers;
layersToPlot = {'ilm' 'nflgcl' 'iplinl' 'inlopl' 'oplonl' 'isos' 'rpe'};% 'rpeSmooth'}; %
% layersToPlot = {'ilm'  'nflgcl'  'inlopl' 'isos' 'rpe'};% 'rpeSmooth'}; %
% layersToPlot = {'ilm' 'isos'};% 'rpeSmooth'}; %

for k = 1:numel(layersToPlot)
    matchedLayers = strcmpi(layersToPlot{k},{retinalLayers(:).name});
    layerToPlotInd = find(matchedLayers == 1);
    retinalLayers1(k)=retinalLayers(layerToPlotInd);
end
retinalLayers=retinalLayers1;
% plot oct image and the obtained retinal layers.
isPlot = 1;
if isPlot,
    imagesc(img);
    axis image; colormap('gray'); hold on; drawnow;
    layersToPlot = {'ilm' 'isos' 'rpe' 'inlopl' 'nflgcl' 'iplinl' 'oplonl'};% 'rpeSmooth'}; %
%     layersToPlot = {'ilm' 'isos' 'rpe' 'inlopl' 'nflgcl'};% 'rpeSmooth'}; %
    hOffset =       [40    0      40    0        0        40       -40      -40]; % for displaying text
    for k = 1:numel(layersToPlot)

        matchedLayers = strcmpi(layersToPlot{k},{retinalLayers1(:).name});
        layerToPlotInd = find(matchedLayers == 1);

        if ~isempty(retinalLayers1(layerToPlotInd).pathX)         
            colora = params.colorarr(k,:);
%             retinalLayers(layerToPlotInd).pathX = smoothdata(retinalLayers(layerToPlotInd).pathX,'gaussian',25);
            plot(retinalLayers1(layerToPlotInd).pathY,retinalLayers1(layerToPlotInd).pathX-1,'-','color',colora,'linewidth',1.5);
            plotInd = round(numel(retinalLayers1(layerToPlotInd).pathX)/2);            
            text(retinalLayers1(layerToPlotInd).pathY(plotInd)+hOffset(k),retinalLayers1(layerToPlotInd).pathX(plotInd)+params.txtOffset,retinalLayers1(layerToPlotInd).name,'color',colora,'linewidth',2);            
            drawnow;
        end % of if ~isempty            

    end % of k
    hold off;
end % of isPlot        
