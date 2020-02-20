

function [rPaths, img,distance] = getRetinalLayersCore(layerName,img,params,rPaths,flag,distance)


if nargin < 3
    display('3 inputs required, getLayers.m');
    return;   
end

szImg = size(img);

switch layerName
    
    case {'roughILMandISOS'}
        
        imgOld = img(:,2:end-1);
        pathsTemp = getHyperReflectiveLayers(imgOld,params.roughILMandISOS);                 
                        
        %save to structure 
        clear rPaths
        rPaths = pathsTemp;
        
        return;        

    case {'ilm' 'isos' 'rpe' 'isos_refine' 'inlopl' 'nflgcl' 'iplinl' 'iplinl_refine' 'oplonl' 'inlopl_refine' 'oplonl_refine'}

        adjMA = params.adjMA;
        adjMB = params.adjMB;
        adjMW = params.adjMW;
        adjMmW = params.adjMmW;        
        
    case {'IF_YOU_WANT_A_SMOOTHER_EDGE_PUT_THE_LAYER_NAMES_HERE'}

        adjMA = params.adjMA;
        adjMB = params.adjMB;
        adjMW = params.adjMWSmo;
        adjMmW = params.adjMmWSmo;
        
end


% initialize region of interest
szImg = size(img);
roiImg = zeros(szImg);

% avoid the top part of image
roiImg(1:20,:) = 0;

% select region of interest based on layers priorly segmented.
        %%ming
        img_show=img;
        startrecord=0;
        distance
        if distance(1)>distance(3)*5 && distance(1)>20
            width_Sunken=150;
        else
            width_Sunken=0;
        end
        width_Sunken
        
for k = 2:szImg(2)-1
    
    switch layerName
        case {'ilm'}
            % define a region (from 'startInd' to 'endInd') near 'ilm'.
            indPathX = find(rPaths(strcmp('ilm',{rPaths.name})).pathY==k);
            
            startInd = rPaths(strcmp('ilm',{rPaths.name})).pathX(indPathX(1)) - params.ilm_0; 
            endInd = rPaths(strcmp('ilm',{rPaths.name})).pathX(indPathX(1)) + params.ilm_1;             
        case {'isos'}            
            % define a region (from 'startInd' to 'endInd') near 'isos'.
            indPathX = find(rPaths(strcmp('isos',{rPaths.name})).pathY==k);            
            
            startInd = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1)) - params.isos_0; 
            endInd = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1)) + params.isos_1;             
            
        case {'rpe'}
            indPathX = find(rPaths(strcmp('isos',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1));
            endInd0 = startInd0+round((rPaths(strcmp('isos',{rPaths.name})).pathXmean-rPaths(strcmp('ilm',{rPaths.name})).pathXmean));

            startInd = startInd0+round(params.rpe_0*(endInd0-startInd0));
            endInd = endInd0-round(params.rpe_1*(endInd0-startInd0));  
            
        case {'isos_refine'}            
            % define a region (from 'startInd' to 'endInd') near 'isos'.
            %%ming
            indPathX = find(rPaths(strcmp('rpe',{rPaths.name})).pathY==k);                        
            startInd0 = rPaths(strcmp('rpe',{rPaths.name})).pathX(indPathX(1));
            endInd0 = rPaths(strcmp('rpe',{rPaths.name})).pathX(indPathX(1));

            startInd = startInd0-40+0*round(params.inlopl_0*(endInd0-startInd0));
            endInd = endInd0-30-0*round(params.inlopl_1*(endInd0-startInd0)); 
        
        case {'inlopl'} %特不准
            % define a region (from 'startInd' to 'endInd') between 'ilm'
            % and 'isos'.
            indPathX = find(rPaths(strcmp('ilm',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('ilm',{rPaths.name})).pathX(indPathX(1));
            indPathX = find(rPaths(strcmp('isos',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1));
                        
            startInd = startInd0+round(params.inlopl_0*(endInd0-startInd0));
            endInd = endInd0-round(params.inlopl_1*(endInd0-startInd0));            
       
        case {'nflgcl'}%准
            % define a region (from 'startInd' to 'endInd') between 'ilm'
            % and 'inlopl'.
            indPathX = find(rPaths(strcmp('ilm',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('ilm',{rPaths.name})).pathX(indPathX(1));            
            
            indPathX = find(rPaths(strcmp('inlopl',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('inlopl',{rPaths.name})).pathX(indPathX(1));

            if width_Sunken~=0
                add= abs(k-distance(2) )/ (szImg(2)-distance(2));
            else
                add=2;
            end
            startInd = startInd0;
            endInd = startInd+params.nflgcl_1*(endInd0-startInd)*add/2;% - round(params.nflgcl_1*(endInd0-startInd0));
  
        case {'iplinl'}
            % define a region (from 'startInd' to 'endInd') between
            % 'nflgcl' and 'inlopl'.
            indPathX = find(rPaths(strcmp('nflgcl',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('nflgcl',{rPaths.name})).pathX(indPathX(1));
            indPathX = find(rPaths(strcmp('inlopl',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('inlopl',{rPaths.name})).pathX(indPathX(1));
%             startInd = startInd0 + round(params.iplinl_0*(endInd0-startInd0));
%             endInd = endInd0 - round(params.iplinl_1*(endInd0-startInd0));
            endInd = endInd0 - 0*round(params.iplinl_1*(endInd0-startInd0));
            add= abs(k-distance(2))/ distance(2); 
            startInd = startInd0 + add*round(params.iplinl_0*(endInd0-startInd0));
            
            if width_Sunken~=0
                if k>=distance(2)-width_Sunken && (k<=distance(2)+width_Sunken)
                    add= abs(k-distance(2))/ (width_Sunken);
                    add=add^(1);
                    endInd = startInd+(endInd0-startInd)*add;% - round(params.nflgcl_1*(endInd0-startInd0));
                end
            end
        case {'iplinl_refine'}
            % define a region (from 'startInd' to 'endInd') between
            % 'nflgcl' and 'inlopl'.
            indPathX = find(rPaths(strcmp('nflgcl',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('nflgcl',{rPaths.name})).pathX(indPathX(1));
            indPathX = find(rPaths(strcmp('inlopl',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('inlopl',{rPaths.name})).pathX(indPathX(1));

            indPathX = find(rPaths(strcmp('iplinl',{rPaths.name})).pathY==k);
            mid0 = rPaths(strcmp('iplinl',{rPaths.name})).pathX(indPathX(1));
            
            startInd=round((startInd0+mid0)/2);
            endInd=round((endInd0+mid0)/2);
            
            if width_Sunken~=0
                if k>=distance(2)-width_Sunken && (k<=distance(2)+width_Sunken)
                    add= abs(k-distance(2))/ (distance(2));
                    endInd = startInd+(endInd0-startInd)*add;% - round(params.nflgcl_1*(endInd0-startInd0));
                end
            end
            
        case {'oplonl'}
            % define a region (from 'startInd' to 'endInd') between
            % 'inlopl' and 'isos'.
            indPathX = find(rPaths(strcmp('inlopl',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('inlopl',{rPaths.name})).pathX(indPathX(1));
            indPathX = find(rPaths(strcmp('isos',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1));

%           startInd = startInd0 + params.oplonl_0;
%           endInd = endInd0 - params.oplonl_1;
            startInd = startInd0 +round(params.oplonl_0*(endInd0-startInd0));
            endInd = endInd0 -round(params.oplonl_1*(endInd0-startInd0));
            
        case {'inlopl_refine'} %
            % define a region (from 'startInd' to 'endInd') between 'ilm'
            % and 'isos'.
            indPathX = find(rPaths(strcmp('iplinl',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('iplinl',{rPaths.name})).pathX(indPathX(1));
            indPathX = find(rPaths(strcmp('oplonl',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('oplonl',{rPaths.name})).pathX(indPathX(1));
            
            endInd = endInd0-0*round(params.inlopl_1*(endInd0-startInd0));
            startInd = startInd0+0*round(params.inlopl_0*(endInd0-startInd0));

            if width_Sunken~=0
                if k>=distance(2)-width_Sunken && (k<=distance(2)+width_Sunken)
                    add= abs(k-distance(2))/ (width_Sunken*2)+0.2;
                    add=add^(1);
        %             endInd = endInd0-0*round(params.inlopl_1*(endInd0-startInd0));
                    endInd = startInd+(endInd0-startInd)*add;% - round(params.nflgcl_1*(endInd0-startInd0));
                end
            end
        case {'oplonl_refine'} %
            % define a region (from 'startInd' to 'endInd') between 'ilm'
            % and 'isos'.
            indPathX = find(rPaths(strcmp('inlopl',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('inlopl',{rPaths.name})).pathX(indPathX(1));
            indPathX = find(rPaths(strcmp('isos',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1));
            
            startInd = startInd0+0*round(params.inlopl_0*(endInd0-startInd0));
            endInd = endInd0-round(0.2*(endInd0-startInd0));
            
            if width_Sunken~=0
                if k>=distance(2)-width_Sunken && (k<=distance(2)+width_Sunken)
                    add= abs(k-distance(2))/ (width_Sunken*2)+0.2;
                    add=add^(1/3);
        %             endInd = endInd0-0*round(params.inlopl_1*(endInd0-startInd0));
                    endInd = startInd+(endInd-startInd)*add;% - round(params.nflgcl_1*(endInd0-startInd0));
                end
            end
    end
    %error checking
    startInd=round(startInd);  
    endInd=round(endInd); 
    if startInd > endInd
        startInd = endInd - 1;
    end       
    if startInd < 1
        startInd = 1;
    end
    if endInd > szImg(1)
        endInd = szImg(1);
    end

               

    % set region of interest at column k from startInd to endInd
    roiImg(startInd:endInd,k) = 1;
%     img_show(startInd:endInd,k) = 1;

    img_show(round(startInd),k) = 1;
    img_show(round(endInd),k) = 1;
    startrecord(k)=startInd;

end

%ensure the 1st and last column is part of the region of interest.

roiImg(:,1)=1;
roiImg(:,end)=1;     

% if (layerName=="inlopl_refine")
%     img_show(:,distance(2))=1;
%     img_show(:,distance(2)-2)=1;
%     img_show(:,distance(2)+2)=1;
% %     imshow(img_show);
% end

if layerName=="isos_refine"
    layerName="isos";
end
if layerName=="inlopl_refine"
    layerName="inlopl";
end
if layerName=="iplinl_refine"
    layerName="iplinl";
end
if layerName=="inlopl_refine"
    layerName="inlopl";
end
if layerName=="oplonl_refine"
    layerName="oplonl";
end

    

        %%寻找lim的凹陷顶点
%         if layerName=="ilm"
%             for iii=2:szImg(2)-1
%                 ddd1=(startrecord(iii)-startrecord(2))/( startrecord(2)-startrecord(szImg(2)-1) );
%                 ddd2=(iii-2)/( 2-(szImg(2)-1) );
%                 ddd3=( 1/(startrecord(2)-startrecord(szImg(2)-1)) )^2;
%                 ddd4=( 1/((2)-(szImg(2)-1)) )^2;
%                 ddd=abs(ddd1-ddd2)/sqrt( ddd3+ddd4 );
%                 distance(3)=distance(3)+ddd;
%                 if distance(1)<ddd
%                     distance(1)=ddd;
%                     distance(2)=iii-5;
%                 end
%             end
%             distance(3)=distance(3)/( szImg(2)-2 );
%         end
        %%寻找lim的凹陷顶点
        if layerName=="ilm"
            d_width=100;
            for iii=2+d_width:szImg(2)-1-d_width
                x1=iii-d_width;
                x2=iii+d_width;
                ddd1=(startrecord(iii)-startrecord(x1))/( startrecord(x1)-startrecord(x2) );

                ddd2=(iii-x1)/( x1-x2 );
                ddd3=( 1/(startrecord(x1)-startrecord(x2)) )^2;
                ddd4=( 1/(x1-x2))^2;
                ddd=abs(ddd1-ddd2)/sqrt( ddd3+ddd4 );
                if startrecord(x1)==startrecord(x2)
                    ddd=startrecord(iii)-startrecord(x1);
                end
                if isnan(ddd)
                    ddd=startrecord(iii)-startrecord(x1);
                end
                distance(3)=distance(3)+ddd;
                if distance(1)<ddd
                    distance(1)=ddd;
                    distance(2)=iii-5;
                end
            end
            distance(3)=distance(3)/( szImg(2)-1-d_width-(2+d_width) );
        end
        


% include only region of interst in the adjacency matrix
includeA = ismember(adjMA, find(roiImg(:) == 1));
includeB = ismember(adjMB, find(roiImg(:) == 1));
keepInd = includeA & includeB;

%     %alternative to ismember, 
%     roiImgOne = find(roiImg(:) == 1)';
%     includeA = sum(bsxfun(@eq,adjMA(:),roiImgOne),2);
%     includeB = sum(bsxfun(@eq,adjMB(:),roiImgOne),2);
%     keepInd = includeA & includeB;

%get the shortestpath
switch layerName
    %bright to dark
    case {'rpe' 'nflgcl' 'oplonl' 'iplinl' }
        adjMatrixW = sparse(adjMA(keepInd),adjMB(keepInd),adjMW(keepInd),numel(img(:)),numel(img(:)));    
        [ ~, path ] = graphshortestpath( adjMatrixW, 1, numel(img(:)) );
        % dist = nan(size(path));
        % for i = 1:numel(path)-1,dist(i)=adjMatrixW(path(i),path(i+1));end
    %dark to bright
    case {'inlopl' 'ilm' 'isos' }
        adjMatrixMW = sparse(adjMA(keepInd),adjMB(keepInd),adjMmW(keepInd),numel(img(:)),numel(img(:)));    
        [ ~, path ] = graphshortestpath( adjMatrixMW, 1, numel(img(:)) );        
        % dist = nan(size(path));
        % for i = 1:numel(path)-1,dist(i)=adjMatrixMW(path(i),path(i+1));end
end

%convert path indices to subscript
[pathX, pathY] = ind2sub(szImg,path);

        %%ming
%         imagesc(img); axis image; colormap('gray'); hold on;
%         plot(pathY,pathX,'r-','linewidth',1); hold on;
%         title(layerName);


%if name layer existed, overwrite it, else add layer info to struct
matchedLayers = strcmpi(layerName,{rPaths(:).name});
layerToPlotInd = find(matchedLayers == 1);
if isempty(layerToPlotInd)
    layerToPlotInd = numel(rPaths)+1;
    rPaths(layerToPlotInd).name = layerName;
end


% save data.
rPaths(layerToPlotInd).path = path;
% rPaths(layerToPlotInd).dist = dist;
rPaths(layerToPlotInd).pathX = pathX;
rPaths(layerToPlotInd).pathY = pathY;
rPaths(layerToPlotInd).pathXmean = mean(rPaths(layerToPlotInd).pathX(gradient(rPaths(layerToPlotInd).pathY)~=0));
%create an additional smoother layer for rpe
isSmoothRpe = 1;
if isSmoothRpe,
%     switch layerName
%         case {'rpe'}  
%             find lines where pathY is on the image
            rpePathInd = gradient(pathY) ~= 0;

            % fit line with cubic smoothing spline
            lambda = 1E-5; %small means really smooth
            pathXpoly = pathX;
            pathYpoly = pathY;
           
%             [pathXpoly(rpePathInd), ~] = csaps(pathY(rpePathInd),pathX(rpePathInd),lambda,pathY(rpePathInd));               
           [pathXpoly(rpePathInd), ~] = csaps(pathY(rpePathInd),pathX(rpePathInd),lambda,pathY(rpePathInd));               
%             [pathXpoly(rpePathInd), ~] = csaps_pt(pathY(rpePathInd),pathX(rpePathInd),...
%                lambda,pathY(rpePathInd));               
            %add layer info to struct
            %layerToPlotInd = numel(rPaths)+1;
            %rPaths(layerToPlotInd).name = 'rpeSmooth';
            %update rpw layer info to struct
            rPaths(layerToPlotInd).pathX = round(pathXpoly);
            rPaths(layerToPlotInd).pathY = round(pathYpoly);            
            rPaths(layerToPlotInd).path = sub2ind(szImg,rPaths(layerToPlotInd).pathX,rPaths(layerToPlotInd).pathY);
            rPaths(layerToPlotInd).pathXmean = mean(rPaths(layerToPlotInd).pathX(gradient(rPaths(layerToPlotInd).pathY)~=0));            

% %         otherwise
% %     end

end