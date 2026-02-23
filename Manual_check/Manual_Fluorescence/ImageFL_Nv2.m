% Version 2.2   @22.11.30   added delSLtr and delCelSL buttons.
% Version 2.1   @22.11.11   added tracing button.

classdef ImageFL_Nv2 < handle    
    properties
        Figure
        Dialog
        Axis
        TransparencyHandle
        Image% this is the background [unit8] sz x 3 image.
        ImageHandle
        MouseIsDown
        PX
        PY
        Scroll
        SList% list of selected spots
        TextX
        TextY
% v1 newly added:
        Map
        LineLayer
        IndLineHidden
        IndSquare
        LineStEnd
        sz
        selSquare
        selLine
        SptCent
        Paths
% v2 newly added:
        Figure2
        Axis2
        ImageHandle2
        Scroll2
        CrossHair
        curSelSpot
        PHMask
        RJ
        CIDinParts  
        HorLineMap  
        HorLineMapOut   
        PHparts     
% v2.1 newly added:
        SptLines    %spot lines
        SptMap      %spot map 
        SLinesMap   %spot lines map
        tempSptMap  %spot line map (for manual correction)
        SptParts    %spot line parts on the tree plot.
        MM          %Mapping method: [4] fast growth, [5] slow growth
        FigSLMap  %table to show SLMap (i.e., SptMap).
        SLmapOVR    %Spot line overide item. For manual correction of SL map.
    end
    
    methods
        function tool = ImageFL_Nv2(J,SList,LineLayerRaw,RJ,indSquare,indLine,LineStEnd,SptCent,pad,PHMaskIn,cell4test,MM,Paths)
            tool.Image=J;%This stores image canvas with RJ, red contours, green x, green spot ID.
            tool.LineLayer=LineLayerRaw;
            tool.IndLineHidden=indLine;
            tool.IndSquare=indSquare;
            tool.LineStEnd=LineStEnd;
            tool.SptCent=SptCent;
            tool.sz=size(J);%static value.
            tool.Paths=Paths;
            tool.selSquare=[];
            tool.selLine=[];
            tool.MouseIsDown = false;
            tool.Map=[];
            tool.CrossHair=pad;
            tool.curSelSpot=[];
            tool.SList=SList;
            tool.PHMask=PHMaskIn;%this is tilized mask [logical] for the current SchN.
            tool.RJ=RJ;
            tool.CIDinParts=cell4test{1,1};%changed on 23.03.03 due to LkgMv2
            tool.PHparts=cell4test{1,2};%changed on 23.03.03 due to LkgMv2
            tool.HorLineMap=cell4test{1,3};
            tool.SptLines=[];
            tool.SptMap=[];
            tool.tempSptMap=[];
            tool.SptParts=[];
            tool.HorLineMapOut=[];
            tool.MM=MM;
            tool.FigSLMap=[];
            tool.SLmapOVR =[];
            tool.SLinesMap=[];
            tool.Figure = figure(...%'MenuBar','none', ...
                                 'NumberTitle','off', ...
                                 'Name','连连看 Pro', ...
                                 'Position',[350, 50, 1500, 300],...% 'CloseRequestFcn',@tool.closeFigure, ...
                                 'WindowButtonMotionFcn', @tool.mouseMove, ...%click down and hold
                                 'WindowButtonDownFcn', @tool.mouseDown, ...%click down
                                 'WindowButtonUpFcn', @tool.mouseUp, ...%release mouse
                                 'WindowKeyPressFcn',@tool.keyPressCB,...
                                 'Resize','on');

            tool.Figure2 = figure(...%'MenuBar','none', ...
                                 'NumberTitle','off', ...
                                 'Name','Raw', ...
                                 'Position',[350, 410, 1500, 300],...
                                 'CloseRequestFcn',@tool.closeFigure, ...%'WindowButtonMotionFcn', @tool.mouseMove, ...%click down and hold
                                 'WindowButtonDownFcn', @tool.mouseDown, ...%click down
                                 'WindowButtonUpFcn', @tool.mouseUp, ...%release mouse
                                 'Resize','on');

            tool.Axis = axes('Parent',tool.Figure,'Position',[0 0 1 1]);
            tool.Axis2 = axes('Parent',tool.Figure2,'Position',[0 0 1 1]);
            %How to generate image for display?
            %out=cat(3,RJ+canvas(:,:,1)+Contour,RJ+2*canvas(:,:,2)+canvas(:,:,1),RJ);%RJ:FLTS; canvas:x-cross and text labels; contour:cellPHMaskContour
            %J=cat(3,RJ+Contour,RJ+2*canvas(:,:,2),RJ);
            %So, out=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
            tool.ImageHandle = imshow(cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3)),'Parent',tool.Axis);%This stores the image (tool.ImageHandle)
            tool.ImageHandle2 = imshow(RJ,'Parent',tool.Axis2);
            tool.Scroll = imscrollpanel(tool.Figure,tool.ImageHandle);
            tool.Scroll2 = imscrollpanel(tool.Figure2,tool.ImageHandle2);

            api = iptgetapi(tool.Scroll);
            api.setMagnificationAndCenter(1,256,100);%<-------------set magnification and center.
            api.addNewLocationCallback(@tool.ScrollFun1);

            api2 = iptgetapi(tool.Scroll2);
            api2.setMagnificationAndCenter(1,256,100);%<-------------set magnification and center.
            api2.addNewLocationCallback(@tool.ScrollFun2);

            loc = api.getVisibleLocation();
            loc=round(loc);
            api2.setVisibleLocation(loc);
            api.setVisibleLocation(loc);


            Jt = ones(tool.sz(1:2));%set up transparent layer.
            Jt = cat(3,Jt,Jt,zeros(size(J,1),size(J,2)));%set up transparent layer.
            hold on
            tool.TransparencyHandle = imshow(Jt,'Parent',tool.Axis2);%This stores the blue layer [0,0,1] (tool.TransparencyHandle)
            hold off
            set(tool.TransparencyHandle, 'AlphaData', zeros(tool.sz(1:2)));
            annotation(tool.Figure,'arrow',[0.0063,0.5],[0.0127,0.6]);

            % initialization of UI components ('Dialog'):
            tool.Dialog =dialog('WindowStyle', 'normal',...
                                'Name', 'ImFL-Master',...
                                'CloseRequestFcn', @tool.closeDialog,...
                                'Position',[10 50 340 300]);
            
            % XY mouse pointer coordinates.
            uipanel(tool.Dialog,'units','pixels','Position',[10,198,130,72]);
            uicontrol(tool.Dialog,'Position',[11,200, 50,12],'Style','text','HorizontalAlignment','right','String','Y:','FontSize',7);
            tool.PY=uicontrol(tool.Dialog,'Position',[62,200, 25,12],'Style','text','HorizontalAlignment','left','String','','FontSize',7);
            uicontrol(tool.Dialog,'Position',[11,213, 50,12],'Style','text','HorizontalAlignment','right','String','X:','FontSize',7);
            tool.PX=uicontrol(tool.Dialog,'Position',[62,213, 25,12],'Style','text','HorizontalAlignment','left','String','','FontSize',7);
%% Press buttons: AddSpot, Conv2DualSpots, DelSpot, GenTree, GenMap, RefreshPic, Skip, SaveList, SavePic, SLsharp, Test(with TableDone),
            % Button Trace
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','TrSpot','Position',[220 105 50 20],'Callback',@tool.btnTraceSpot);
            % Button delSLtrace
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','delSLtr','Position',[220 155 50 20],'Callback',@tool.btnRemSLtrace);
            % Button delCellSLines
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','delCelSL','Position',[220 130 50 20],'Callback',@tool.btnRemCellSL);




%             % Button AddSpot:
%             uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','AddSpot','Position',[280 105 50 20],'Callback',@tool.btnAddSpot);
            tool.TextX = uicontrol('Parent',tool.Dialog,'Style','edit','Units','pixels','Position',[210 105 25 20],'Visible','off');
            tool.TextY = uicontrol('Parent',tool.Dialog,'Style','edit','Units','pixels','Position',[250 105 25 20],'Visible','off');
%             uicontrol('Parent',tool.Dialog,'Style','text','Units','pixels','String','X','Position',[200 103 10 20]);
%             uicontrol('Parent',tool.Dialog,'Style','text','Units','pixels','String','Y','Position',[240 103 10 20]);
            % Button 1->3: Convert 2 dual spots:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','1->3','Position',[280 155 50 20],'Callback',@tool.btnConv2DualSpots);
            % Button Delete spot:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','DelSpot','Position',[280 130 50 20],'Callback',@tool.btnDelSpot);
            % Button GenTree: generate a tree and see if everything is OK
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','genTree','Position',[220 180 50 40],'Callback',@tool.btnGenTree);
            % Button GenMap:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','genMap','Position',[280 55 50 20],'Callback',@tool.btnGenMap);
            % Button Refresh background pictures:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','Refresh','Position',[280,220,50,40],'Callback',@tool.btnRefreshPic);
            % Button Skip: skip this side channel
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','Skip','Position',[220 10 50 40],'Callback',@tool.btnSkip);  
            % Button SaveList button:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','SaveList','Position',[220 80 50 20],'Callback',@tool.btnSaveList);
            % Button SavePic button:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','SavePic','Position',[280 80 50 20],'Callback',@tool.btnSavePic);
            % SL#: show the spot line index at the initial spot.
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','SL#','Position',[220 55 50 20],'Callback',@tool.btnSLsharp);
            % Button test: test if every line is correctly linked:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','Test','Position',[280 180 50 40],'Callback',@tool.btnTest);
            % Button Done:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','Done','Position',[280 10 50 40],'Callback',@tool.btnDone);

%             tool.btnTest;

            uiwait(tool.Dialog)

            


        end
%% GUI functions:
% Press buttons: AddSpot, Conv2DualSpots, DelSpot, GenTree, GenMap, RefreshPic, Skip, SaveList, SavePic, SLsharp, Test(with TableDone),
        function btnTraceSpot(tool,~,~) %trace-draw the line in this div cyc. [q]
            %1. select a spot to be traced.
            tool.btnGenMap;
            if isempty(tool.curSelSpot)
                disp('Select a spot first!');
                return
            end
            curSelSpt=tool.curSelSpot;

            %2. retrieve the PHpart containing this spot
            curSelCID=tool.SList(curSelSpt,1:2);
            selPHId=cellfun(@(x) ~isempty(find(x(:,1)==curSelCID(1) & x(:,2)==curSelCID(2), 1)),tool.CIDinParts);
            stIdinPHpart=find(tool.CIDinParts{selPHId}(:,1)==curSelCID(1));
            remCID=tool.CIDinParts{selPHId}(stIdinPHpart:end,:);

            %3. gather sptIDs that appear in each remaining CID
            sptIdInEachCID=cell(size(remCID,1),1);
            for i=1:size(remCID,1)
                sptIdInEachCID{i}=find(tool.SList(:,1)==remCID(i,1) & tool.SList(:,2)==remCID(i,2));
            end
            
            %4. if the selected spot is the only one, link the most likely spots in the PHpart
            if isscalar(sptIdInEachCID{1})
                route=ImgFLN_traceSpotInRemCID(sptIdInEachCID,tool.SList,curSelSpt,tool.Map,[1,0]);
            elseif length(sptIdInEachCID{1})==2
                route=ImgFLN_traceSpotInRemCID(sptIdInEachCID,tool.SList,curSelSpt,tool.Map,[2,0]);
%             else
            end
            %5. draw new lines.
            if isempty(route)
                disp(['No new lines from #',num2str(curSelSpt),' can be projected.']);
                return;
            end
            disp(['Projecting new lines from #',num2str(curSelSpt),'.']);
            [tool.LineLayer,newLineIdx]=ImgFLN_createNewLineBatch(tool.LineLayer,route,tool.SptCent);
            tool.IndLineHidden=[tool.IndLineHidden;newLineIdx];
            newLineStEnd=[route,tool.SptCent(route(:,1),:),tool.SptCent(route(:,2),:)];
            tool.LineStEnd=[tool.LineStEnd;newLineStEnd];
            out=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
            set(tool.ImageHandle,'CData',out);
            tool.selSquare=[];
            tool.curSelSpot=[];
        end

        function btnAddSpot(tool,~,~)% adding a new spot
            tool.btnGenMap;
            X=str2double(tool.TextX.String);
            Y=str2double(tool.TextY.String);
            newSptName=length(tool.IndSquare)+1;
            
            tool.IndSquare=[tool.IndSquare;{Cent2Idx([X,Y],tool.sz(1:2))}];
            tool.SptCent=[tool.SptCent;[X,Y]];
            tool.Map=[tool.Map;{[]}];
            out=tool.Image;
            rawLength=size(tool.SList,1);
            [out,tool.SList]=drawSptMapFLvImgFL_Nintv2(out,tool.SptCent,tool.Map,tool.SList,tool.PHMask,1);
            newLength=size(tool.SList,1);
            if rawLength==newLength
                tool.IndSquare(end,:)=[];
                tool.SptCent(end,:)=[];
                tool.Map(end,:)=[];
                newSptName=0;
            end
            tool.Image=out;
            J=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
            set(tool.ImageHandle,'CData',J);
            disp(['Added Spot #',num2str(newSptName)]);
        end  %done

        function btnConv2DualSpots(tool,~,~)%1->3 convert 1 spot to 3.
            if isempty(tool.Map)
                tool.btnGenMap;
            end
            if isempty(tool.curSelSpot)
                disp('Select a spot first!');
                return
            end
            curSelSpt=tool.curSelSpot;
            tool.curSelSpot=[];

            %add two new spots
            X=tool.SptCent(curSelSpt,1);
            Y=tool.SptCent(curSelSpt,2);
            newSptNames=[length(tool.IndSquare)+1,length(tool.IndSquare)+2];%keep the original one and add two additional ones [up,low].
            tool.IndSquare=[tool.IndSquare;{Cent2Idx([X,Y-4],tool.sz(1:2))};{Cent2Idx([X,Y+4],tool.sz(1:2))}];
            tool.SptCent=[tool.SptCent;[X,Y-4];[X,Y+4]];
            tool.Map=[tool.Map;{[]};{[]}];

            %remove the mapping info for the old spot:
            tool.Map(curSelSpt,:)={[]};

            %remove the line (if any) that starts from the old spot:
            del=find(tool.LineStEnd(:,1)==curSelSpt);
            if ~isempty(del)
                for i=1:length(del)
                    disp(['Deleted line between #',num2str(tool.LineStEnd(del(i),1)),' to #', num2str(tool.LineStEnd(del(i),2)),'.']);
                    [tool.LineLayer]=ImgFLN_removeOldLine(tool.LineLayer,tool.LineStEnd(del(i),:));


                end
                tool.IndLineHidden(del)=[];
                tool.LineStEnd(del,:)=[];
                tool.selLine=[];
            end
            out=tool.Image;
            [out,tool.SList]=drawSptMapFLvImgFL_Nintv2(out,tool.SptCent,tool.Map,tool.SList,tool.PHMask,2);
            tool.Image=out;
            J=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
            set(tool.ImageHandle,'CData',J);
            disp(['Terminated a Spot Line at #',num2str(curSelSpt)]);
            disp(['Added two new spots #',num2str(newSptNames(1)),' and #',num2str(newSptNames(2))]);
        end

        function btnDelSpot(tool,~,~)%deleting an existing spot
    %1. remove entry in SList
            if isempty(tool.Map)
                tool.btnGenMap;
            end
            if isempty(tool.curSelSpot)
                disp('Select a spot first!');
                return
            end
            curSelSpt=tool.curSelSpot;
            tool.SptCent(curSelSpt,:)=nan(1,2);
            tool.curSelSpot=[];
            tool.SList(curSelSpt,:)=nan(1,9);

    %2. remove its own entry in Map, as well as any line that points to this spot.
            tool.Map(curSelSpt)={[]};
            idx=find(cell2mat(cellfun(@(x) ismember(curSelSpt,x),tool.Map,'UniformOutput',false)));
            for i1=1:length(idx)
                vx=tool.Map{idx(i1),1};
                vx(vx==curSelSpt)=[];
                tool.Map{idx(i1),1}=vx;
            end

    %3. remove the indSquare
            tool.IndSquare(curSelSpt,1)={[]};
            
    %4. remove the line from the image.
            del=find(tool.LineStEnd(:,1)==curSelSpt | tool.LineStEnd(:,2)==curSelSpt);
            if ~isempty(del)
                for i=1:length(del)
                    disp(['Deleted line between #',num2str(tool.LineStEnd(del(i),1)),' to #', num2str(tool.LineStEnd(del(i),2)),'.']);
                    [tool.LineLayer]=ImgFLN_removeOldLine(tool.LineLayer,tool.LineStEnd(del(i),:));
                end
                tool.IndLineHidden(del)=[];
                tool.LineStEnd(del,:)=[];
                tool.selLine=[];
            end

  %5. remove from the image:
            tool.Image=ImgFLN_removeOldSpot(tool.RJ,tool.SptCent,tool.PHMask);
            J=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
            set(tool.ImageHandle,'CData',J);
            disp(['Deleted Spot #',num2str(curSelSpt)]);
  %6. (Optional) validate if the deletion is thorough.
  
        end

        function btnGenTree(tool,~,~)% generating a tree with both blue and orange lines.
            if isempty(tool.SLinesMap)
                clc;
                disp('Fix problems listed in Test.');
                return;
            end
            ImgFLN_drawSPTtree(tool.SptParts,tool.SLinesMap,0.05,'visible',[],tool.Paths);
        end
        
        function btnGenMap(tool,~,~)% generating a spot map from spots and connections. Called in other functions.
            map=cell(length(tool.IndSquare),1);
            for i=1:size(tool.LineStEnd,1)
                A=map{tool.LineStEnd(i,1),1};
                A=cat(2,A,tool.LineStEnd(i,2));
                map{tool.LineStEnd(i,1),1}=A;
            end
            disp('Map created.');
            tool.Map=map;
        end  %done

        function btnRefreshPic(tool,~,~) %Refresh the background image by re-drawing the line layer.
            canvas=uint8(zeros(tool.sz(1:2)));
            for iLin=1:size(tool.LineStEnd,1)
                stPoint=tool.LineStEnd(iLin,3:4);
                endPoint=tool.LineStEnd(iLin,5:6);
                canvas=insertShape(canvas,'Line',[stPoint,endPoint],'LineWidth',1,'Color','red');
                canvas(stPoint(1,2)-4:stPoint(1,2)+4,stPoint(1,1)-4:stPoint(1,1)+4,1)=0;
                canvas(endPoint(1,2)-4:endPoint(1,2)+4,endPoint(1,1)-4:endPoint(1,1)+4,1)=0;
            end
            tool.LineLayer=canvas(:,:,1);
            out=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
            set(tool.ImageHandle,'CData',out);
            disp('Line map refreshed.')
        end  %done

        function btnSaveList(tool, ~,~)%temperarily save the SList and spot map.
            tool.btnGenMap;
            MapX=tool.Map;
            SListX=tool.SList;
            SptCentX=tool.SptCent;
            L=cellfun(@(x) length(x),MapX);
            MapN=nan(length(L),3);
            for i=1:size(MapX,1)
                if L(i)~=0
                    MapN(i,1:L(i))=MapX{i,1};
                end
            end
            load(tool.Paths{2,1});
            SchN=tool.Paths{3,1};

            SptMap{1,SchN}=MapN;
            SptList{1,SchN}=SListX;
            SptCent{1,SchN}=SptCentX;
            save(tool.Paths{2,1},'SptMap','SptList','SptCent','-append');
            disp('Spot Map List saved.')
        end
        
        function btnSavePic(tool,~,~)%save the picture with spot lines, spots and cell contours.
            out=get(tool.ImageHandle,'CData');
            imwrite(out,tool.Paths{1,1});
        end  %done

        function btnSkip(tool, ~,~)%Skip this side channel
            tool.SList=[];
            tool.SptCent=[];
            tool.CIDinParts=[];
            tool.HorLineMap=[];
            tool.PHparts=[];
            tool.SptLines=[];
            tool.SptMap=[];
            tool.SptParts=[];
            delete(tool.Figure2);
            delete(tool.Dialog);%<-----perserve permanently
            delete(tool.Figure);%<-----perserve permanently   
        end

        function btnSLsharp(tool,~,~)%show SL index, manually correct SL map and test.
        %1. Add SL number in front of each initial spot.
            tool.btnGenMap;
            curMap = tool.reshapeMap();
            sptPartsInfo=ImgFLN_sptMap2sptParts_core(curMap,tool.SList);
            SLines=ImgFLN_MarkSL(sptPartsInfo,tool.CIDinParts);
            Q=cellfun(@(x) x(1,1), SLines);
            locList=tool.SptCent(Q,:)-[14,0];
            canvas=tool.ImageHandle.CData;
            for iLin=1:length(Q)
                if locList(iLin,1)<0
                    locList(iLin,1)=0;
                end
                canvas = insertText(canvas,locList(iLin,:), num2str(iLin),'BoxOpacity',0,'FontSize',12,'TextColor','yellow');
            end
            set(tool.ImageHandle,'CData',canvas);

        %2. display tempSLMap for manual correction: (call up a GUI with a test button)
            tool.FigSLMap = figure('Position',[10 390 340 720],'Name','Spot Line Map');
            if isempty(tool.tempSptMap)
                dat=zeros(1,3);
            else
                dat=[(1:1:size(tool.tempSptMap,1))',tool.tempSptMap];
            end
            cnames = {'Mom','DaugtherUp','DaughterLow'};
            t = uitable('Parent',tool.FigSLMap,'Data',dat,'ColumnName',cnames,'Position',[10 10 280 710]);                        
            set(t,'ColumnEditable',true);
            uicontrol('Parent',tool.FigSLMap,'Style','pushbutton','String','Done','Position',[292 10 38 30],'Callback',@tool.btnTableDone);
            uiwait(tool.FigSLMap);
            %3. redo test    
            curMap = tool.reshapeMap();
            sptPartsInfo=ImgFLN_sptMap2sptParts_core(curMap,tool.SList);
            [curSptParts,curSptMap,~,curSptLines,curHorLineMap]=ImgFLN_testSptDrawing(tool.tempSptMap,sptPartsInfo,tool.HorLineMap,tool.CIDinParts,tool.PHparts,tool.MM);
            if ~isempty(curSptParts)
                tool.SptParts=curSptParts;
                tool.SptMap=curMap;
                tool.SptLines=curSptLines;
                tool.SLinesMap=curSptMap;
                tool.HorLineMapOut=curHorLineMap;
                tool.btnSaveList;
            end
            if 1%meaning the map is problematic
                
            end
        end

        function btnRemSLtrace(tool,~,~)% remove all lines connecting to the selected spot.[f]
            tool.btnGenMap;
            if isempty(tool.curSelSpot)
                disp('Select a spot first!');
                return
            end
            curSelSpt=tool.curSelSpot;
            curMap = tool.reshapeMap();
            sptPartsInfo=ImgFLN_sptMap2sptParts_core(curMap,tool.SList);%this generate Slines.

            remSptIndList=ImgFLN_genremSptIndList(curSelSpt,curMap);%remove spot line from current selected point.
            del=[];
            for i=1:size(remSptIndList,1)
                tool.selLine=find(tool.LineStEnd(:,1)==remSptIndList(i,1) & tool.LineStEnd(:,2)==remSptIndList(i,2),1);
                [tool.LineLayer]=ImgFLN_removeOldLine(tool.LineLayer,tool.LineStEnd(tool.selLine,:));
                del=cat(1,del,tool.selLine);
            end
            %remove this line from hidden layer
            tool.IndLineHidden(del,:)=[];
            tool.LineStEnd(del,:)=[];
            tool.selLine=[];
            out=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
            set(tool.ImageHandle,'CData',out);
        end

        function btnRemCellSL(tool,~,~) % remove all lines from current div cycle after the selected time point. [d]
            if isempty(tool.Map)
                tool.btnGenMap;
            end
            if isempty(tool.curSelSpot)
                disp('Select a spot first!');
                return
            end
            
            curSelSpt=tool.curSelSpot;
            curSelCID=tool.SList(curSelSpt,1:2);%which cell is current selected spot in?
            curSelPHLineIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==curSelCID(1) & x(:,2)==curSelCID(2),1)),tool.CIDinParts),1);
            curSelPHLCIDs=tool.CIDinParts{curSelPHLineIdx}(tool.CIDinParts{curSelPHLineIdx}(:,1)>=curSelCID(1),:);
            remSptIndList=[];
            for i=1:size(curSelPHLCIDs,1)%for each CID in this PHline
                selSptIdx=find(tool.SList(:,1)==curSelPHLCIDs(i,1) & tool.SList(:,2)==curSelPHLCIDs(i,2));
                if ~isempty(selSptIdx)

                    for j=1:length(selSptIdx)
                        remSptIndList=[remSptIndList;tool.LineStEnd(tool.LineStEnd(:,1)==selSptIdx(j),1:2)];
                    end
                end
            end
            del=[];
            for i=1:size(remSptIndList,1)
                tool.selLine=find(tool.LineStEnd(:,1)==remSptIndList(i,1) & tool.LineStEnd(:,2)==remSptIndList(i,2),1);
                [tool.LineLayer]=ImgFLN_removeOldLine(tool.LineLayer,tool.LineStEnd(tool.selLine,:));
                del=cat(1,del,tool.selLine);
            end
            %remove this line from hidden layer
            tool.IndLineHidden(del,:)=[];
            tool.LineStEnd(del,:)=[];
            out=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
            set(tool.ImageHandle,'CData',out);
            tool.selLine=[];

        end

        function btnTest(tool,~,~)% test if the current connection is valid.
            tool.btnGenMap;
            curMap = tool.reshapeMap();
            sptPartsInfo=ImgFLN_sptMap2sptParts_core(curMap,tool.SList);
            [curSptParts,curSLinesMap,curTempSptMap,curSptLines,curHorLineMap,errSpots,sw]=ImgFLN_testSptDrawing([],sptPartsInfo,tool.HorLineMap,tool.CIDinParts,tool.PHparts,tool.MM);
            %if curSptLines has problems, sw is false curSptParts is empty. In this case, provide an option to manually correct the SLine Map and test.
            if ~isempty(curSptParts)%all good. proceed to done.
                tool.SptParts=curSptParts;
%                 tool.SptMap=curSptMap;
                tool.SptMap=curMap;
                tool.SptLines=curSptLines;
                tool.SLinesMap=curSLinesMap;
                tool.HorLineMapOut=curHorLineMap;
                tool.btnSaveList;
            else
                tool.SLinesMap=[];
                %V2.03: mark wrong spot lines with red marker. May just draw red circles on the axis. When clicking inside the circle, remove it.
%                 set(tool.ImageHandle,'CData',out);
            end
            if ~sw%if the spot line map is problematic, sw is false. This part will pop-up a window with SLmap for manual correction.
                tool.tempSptMap=curTempSptMap;
            end
        end

        function btnTableDone(tool,~,~)% attach to the btnSLsharp.
            dat=get(tool.FigSLMap.Children(2,1),'Data');
            if ~isempty(tool.tempSptMap)
                tool.tempSptMap=dat(:,2:3);
                tool.SLinesMap=dat(:,2:3);
            else
                tool.SLmapOVR=[tool.SLmapOVR;dat];
            end
            delete(tool.FigSLMap);
            disp('Spot line map updated.')
        end
        
        function keyPressCB(tool,~,event)% place to go if want to set key shortcuts.
            switch event.Key
                case 'a'%transfer current mouse coordinates to the edit box.
                    
                    set(tool.TextX,'String',tool.PX.String);
                    set(tool.TextY,'String',tool.PY.String);
                    tool.btnAddSpot;
                case 'q'
                    tool.btnTraceSpot;
                case 'd'
                    tool.btnRemCellSL
                case 'f'
                    tool.btnRemSLtrace
            end
        end  %done

% original housekeeping functions:
        function mouseDown(tool,~,~)
            %search in indSquare and indLine. If more than one selected,
            %see which one is closer to the center.
            tool.MouseIsDown = true;
            p = tool.Axis.CurrentPoint;
            if p(1,2)<=0 || p(1,1)<=0 || p(1,2)>=tool.sz(1) || p(1,1)>=tool.sz(2)
                return;
            end
            mouseInd=sub2ind(tool.sz(1:2),round(p(1,2)),round(p(1,1)));
            %check which button was pressed: From left to right: 'normal','extend','alt'.
            tool.selSquare=find(cellfun(@(x) ismember(mouseInd,x),tool.IndSquare));
            tool.selLine=find(cellfun(@(x) ismember(mouseInd,x),tool.IndLineHidden));
            if length(tool.selSquare)>1
                tool.selSquare=tool.selSquare(1);
            end
            if length(tool.selLine)>1
                tool.selLine=[];
                disp('selected two or more lines.')
            end

        end% done

        function mouseUp(tool,~,~)
            tool.MouseIsDown = false;
        %if only left button is released, do these：
            kkk=get(tool.Figure,'SelectionType');
            if strcmp(kkk,'normal')
                p = tool.Axis.CurrentPoint;
                if p(1,2)<=0 || p(1,1)<=0 || p(1,2)>=tool.sz(1) || p(1,1)>=tool.sz(2)
                    return;
                end
                mouseInd=sub2ind(tool.sz(1:2),round(p(1,2)),round(p(1,1)));
                releaseSelSquare=find(cellfun(@(x) ismember(mouseInd,x),tool.IndSquare));
                if isempty(releaseSelSquare) || isempty(tool.selSquare)
                    return;
                end
                if releaseSelSquare~=tool.selSquare%different selSquare when Up from when Down
                    disp(['Created new line between #',num2str(tool.selSquare),' to #', num2str(releaseSelSquare),'.']);
                    [tool.LineLayer,newLineIdx]=ImgFLN_createNewLine(tool.LineLayer,[tool.SptCent(tool.selSquare,:),tool.SptCent(releaseSelSquare,:)]);
                    tool.IndLineHidden=[tool.IndLineHidden;newLineIdx];
                    if tool.SptCent(tool.selSquare,1)<tool.SptCent(releaseSelSquare,1)%always from earlier spots to later spots
                        tool.LineStEnd=[tool.LineStEnd;[tool.selSquare,releaseSelSquare,tool.SptCent(tool.selSquare,:),tool.SptCent(releaseSelSquare,:)]];
                    else
                        tool.LineStEnd=[tool.LineStEnd;[releaseSelSquare,tool.selSquare,tool.SptCent(releaseSelSquare,:),tool.SptCent(tool.selSquare,:)]];
                    end
                    out=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
                    set(tool.ImageHandle,'CData',out);
                    tool.selSquare=[];
                else
                    tool.curSelSpot=tool.selSquare;
                    disp(['Selected #',num2str(tool.selSquare),' but did nothing.']);
                    tool.selSquare=[];
                end
            end
        %if only right button is released, do these：     
            if strcmp(kkk,'alt')
                if ~isempty(tool.selLine)
                    disp(['Deleted line between #',num2str(tool.LineStEnd(tool.selLine,1)),' to #', num2str(tool.LineStEnd(tool.selLine,2)),'.']);
                    [tool.LineLayer]=ImgFLN_removeOldLine(tool.LineLayer,tool.LineStEnd(tool.selLine,:));

                    %remove this line from hidden layer
                    tool.IndLineHidden(tool.selLine)=[];
                    tool.LineStEnd(tool.selLine,:)=[];
                    out=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
                    set(tool.ImageHandle,'CData',out);
                    tool.selLine=[];
                end
            end
        end  %done
        
        function mouseMove(tool,~,~)
            if tool.MouseIsDown
                return;
            end
            p = tool.Axis.CurrentPoint;%<-current coor of pointer.
            col = round(p(1,1));
            row = round(p(1,2));
            set(tool.PX,'String',num2str(col));
            set(tool.PY,'String',num2str(row));
            if row>4 && col>4 && row<tool.sz(1)-4 && col<tool.sz(2)-4 && ~tool.MouseIsDown
                fullMask=false(tool.sz(1:2));
                fullMask((row-4):(row+4),(col-4):(col+4))=tool.CrossHair;
                set(tool.TransparencyHandle, 'AlphaData', fullMask);
            end
        end

        function ScrollFun1(tool,~)
                api1 = iptgetapi(tool.Scroll);
                loc = api1.getVisibleLocation();
                loc=round(loc);
                api2 = iptgetapi(tool.Scroll2);
                api2.setVisibleLocation(loc);
        end
        function ScrollFun2(tool,~)
                api2 = iptgetapi(tool.Scroll2);
                loc = api2.getVisibleLocation();
                loc=round(loc);
                api1 = iptgetapi(tool.Scroll);
                api1.setVisibleLocation(loc);
        end
        
        function closeDialog(tool,src,callbackdata)%<----keep
            delete(tool.Figure2);
            delete(tool.Dialog);
            delete(tool.Figure);
        end
        
        function btnDone(tool,src,callbackdata)%<----working on
            tool.btnGenMap;
            delete(tool.Figure2);
            delete(tool.Dialog);%<-----perserve permanently
            delete(tool.Figure);%<-----perserve permanently
        end
        
        function closeFigure(tool,src,callbackdata)%<----keep
            delete(tool.Figure2);
            delete(tool.Figure);
            delete(tool.Dialog);
        end
%% nested functions:

        function [out]=reshapeMap(in)
            in=in.Map;        
            out=nan(length(in),3);
            for i=1:length(in)
                curIn=unique(in{i,1});
                A=length(curIn);
                if A~=0 && A<=3
                    out(i,1:A)=curIn;
                elseif A~=0 && A>3
                    out(i,1:3)=curIn(1:3);
                end
            end
        end
        function [uit]=getSptIdxInEachCID(in,xy)
            CIDs=in.CIDinParts{selPHId}(stIdinPHpart:end,:);
            uit=cell(size(CIDs,1),1);
            for i=1:size(CIDs,1)
                uit{i}=find(SList(:,1)==CIDs(i,1) & SList(:,1)==CIDs(i,2));
            end
        end
        

    end
end
