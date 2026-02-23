classdef imgABv3 < handle  
    %Version v3.0   @25.06.09   added mutual bottom function.
    %Version v2.1   @23.04.20   Fixed an issue during mesh fitting.
    %Version v2.0   @23.03.06   Output Infcellv2 containing contour fitting data.
    %                           Now, Ulist1 stores lineage problem while Ulist2 stores fitting problem.
    %                           The test is passed when all cells in selected lineage is well fitted.
    properties
        Axis
        DelAfter
        DelUnder
        Dialog
        DidAnnotate
        Figure
        isWasInUlist %if it is freshly entering fitting. Used in btnTest.
        Image
        ImageHandle
        TransparencyHandle
        LabelIndex
        LabelMasks      %mask right before Button Test. Manually modified.
        LabelMasksPrev  %mask after Button Test. For generating modification.
        LineageToSave
        LowerThreshold
        minAL           %minimum area and length, defined in ManCheckStageII.
        meshTempStore
        MouseIsDown
        NLabels
        PenSize
        PenSizeText
        prevICv2
        prevUlist2
        prevUmsk2
        PX %Pointer coordinate X
        PY %Pointer coordinate Y
        RadioDraw
        RadioErase
        RawMsk          %mask from 'predict'. Used for Button Reset
        Scroll
        Slider
        SliderMin
        SliderMax
        UpperThreshold
        U1 %Ulist row number

    end
    
    methods
        
        function tool = imgABv3(I,nLabels,RawMsk,TxtIn,minArea,minLength,minY,cutEndY,varargin)
            tool.LineageToSave={};
            tool.DidAnnotate = 0;
            tool.LabelIndex = 1;
            tool.Image = I;
            tool.minAL=[minArea,minLength,minY];
            J = ones(size(tool.Image));%set up transparent layer.
            J = cat(3,zeros(size(I,1),size(I,2),2),J);%set up transparent layer.
            tool.NLabels = nLabels;
            tool.Figure = figure(...%'MenuBar','none', ...
                                 'NumberTitle','off', ...
                                 'Name','Image', ...
                                 'Position',[350, 50, 1500, 768],...
                                 'CloseRequestFcn',@tool.closeFigure, ...
                                 'WindowButtonMotionFcn', @tool.mouseMove, ...%click down and hold
                                 'WindowButtonDownFcn', @tool.mouseDown, ...%click down
                                 'WindowButtonUpFcn', @tool.mouseUp, ...%release mouse
                                 'WindowKeyPressFcn',@tool.keyPressCB,...
                                 'Resize','on');

            tool.Axis = axes('Parent',tool.Figure,'Position',[0 0 1 1]);
            tool.ImageHandle = imshow(tool.Image);%This stores the image (tool.ImageHandle)
            text(tool.Axis,10,10,TxtIn{1,1},'FontSize',20,'Color','y');%displays side-channel number
            text(tool.Axis,TxtIn{1,2},TxtIn{1,3},TxtIn{1,4},'FontSize',12,'Color','g');%displays time point number
            tool.Scroll = imscrollpanel(tool.Figure,tool.ImageHandle);
            api = iptgetapi(tool.Scroll);
            api.setMagnificationAndCenter(3,256,100);%<-------------set magnification and center.
            hold on
            tool.TransparencyHandle = imshow(J);%This stores the blue layer [0,0,1] (tool.TransparencyHandle)
            hold off
            tool.LabelMasksPrev=RawMsk;%for 1st time enter.
            tool.RawMsk=RawMsk;
            tool.LabelMasks = RawMsk;
            tool.LabelIndex = 1;
            set(tool.TransparencyHandle, 'AlphaData', zeros(size(tool.Image)));%<------modify here? Initially, set transparency to 0.
            tool.MouseIsDown = false;
            tool.Dialog = dialog('WindowStyle', 'normal',...
                                'Name', 'IAB',...
                                'CloseRequestFcn', @tool.closeDialog,...
                                'Position',[50 50 300 300],...
                                'WindowKeyPressFcn',@tool.keyPressCB);
            labels = cell(1,nLabels);
            for i = 1:nLabels
                labels{i} = sprintf('Class %d',i);
            end

            tool.isWasInUlist=true;
            tool.meshTempStore=[];
            tool.prevUlist2=[];
            tool.prevUmsk2=[];
            tool.prevICv2=[];

            %Static text panel
            uipanel(tool.Dialog,'units','pixels','Position',[10,198,130,72]);
            uicontrol(tool.Dialog,'Position',[11,200, 50,12],'Style','text','HorizontalAlignment','right','String','Y:','FontSize',7);
            tool.PY=uicontrol(tool.Dialog,'Position',[62,200, 25,12],'Style','text','HorizontalAlignment','left','String','','FontSize',7);
            uicontrol(tool.Dialog,'Position',[11,213, 50,12],'Style','text','HorizontalAlignment','right','String','X:','FontSize',7);
            tool.PX=uicontrol(tool.Dialog,'Position',[62,213, 25,12],'Style','text','HorizontalAlignment','left','String','','FontSize',7);
            
            uicontrol(tool.Dialog,'Position',[11,226, 50,12],'Style','text','HorizontalAlignment','right','String','U-List#:','FontSize',7);
            uicontrol(tool.Dialog,'Position',[11,240, 50,12],'Style','text','HorizontalAlignment','right','String','Pos:','FontSize',7);
            uicontrol(tool.Dialog,'Position',[62,240, 20,12],'Style','text','HorizontalAlignment','left','String',TxtIn{1,5},'FontSize',7);
            uicontrol(tool.Dialog,'Position',[11,254, 50,12],'Style','text','HorizontalAlignment','right','String','SCh#:','FontSize',7);
            uicontrol(tool.Dialog,'Position',[62,254, 20,12],'Style','text','HorizontalAlignment','left','String',TxtIn{1,6},'FontSize',7);
            tool.U1=uicontrol(tool.Dialog,'Position',[62,226, 40,12],'Style','text','HorizontalAlignment','left','String','','FontSize',7);
            uicontrol(tool.Dialog,'Position',[10,271, 100,12],'Style','text','HorizontalAlignment','left','String','v3.0 ZZ@260609','FontSize',7,'ForegroundColor',[0.5,0.5,0.5]);

            % TEST button:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','TEST [t]','Position',[160 230 130 40],'Callback',@tool.btnTESTPush);

            % Mutual Bottom button
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','Bottom [b]','Position',[160 190 60 20],'Callback',@tool.btnBotPush);
            
            % Deterministic controls:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','Emptize','Position',[160 120 130 20],'Callback',@tool.btnEmptizePush);
            tool.DelAfter=uicontrol('Parent',tool.Dialog,'Style','edit','String','1','Position',[230 170 60 20]);
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','DelAfter','Position',[160 170 60 20],'Callback',@tool.btnDelAfter);
            tool.DelUnder=uicontrol('Parent',tool.Dialog,'Style','edit','String',num2str(cutEndY),'Position',[230 150 60 20]);
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','DelUnder','Position',[160 150 60 20],'Callback',@tool.btnDelUnder);          

            % Erase/Draw
            tool.RadioDraw = uicontrol('Parent',tool.Dialog,'Style','radiobutton','Position',[10 170 70 20],'String','Draw [s]','Callback',@tool.radioDraw);
            tool.RadioErase = uicontrol('Parent',tool.Dialog,'Style','radiobutton','Position',[10 150 70 20],'String','Erase [s]','Callback',@tool.radioErase);
            tool.RadioDraw.Value = 1;
            
            % pencil/eraser slider
            uicontrol('Parent',tool.Dialog,'Style','text','String','Pencil/Eraser Size','Position',[10 120 100 20],'HorizontalAlignment','left');
            tool.PenSizeText = uicontrol('Parent',tool.Dialog,'Style','text','String','5','Position',[125 90 50 20],'HorizontalAlignment','center');
            uicontrol('Parent',tool.Dialog,'Style','edit','String','10','Position',[240 90 50 20],'HorizontalAlignment','right','Callback',@tool.changeSliderRange,'Tag','sliderMax');
            uicontrol('Parent',tool.Dialog,'Style','edit','String','1','Position',[10 90 50 20],'HorizontalAlignment','left','Callback',@tool.changeSliderRange,'Tag','sliderMin');
            tool.PenSize = 4;
            tool.Slider = uicontrol('Parent',tool.Dialog,'Style','slider','Min',1,'Max',10,'Value',tool.PenSize,'Position',[10 70 280 20],'Callback',@tool.sliderManage,'Tag','pss');
            addlistener(tool.Slider,'Value','PostSet',@tool.continuousSliderManage);

            % done button
            if isempty(varargin)
                doneButtonLabel = 'Done';
            else
                doneButtonLabel = varargin;
            end
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String',doneButtonLabel,'Position',[10 10 130 40],'Callback',@tool.buttonDonePushed);
            
            % Restart button
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','Restart','Position',[160 10 130 40],'Callback',@tool.btnRstPushed);
            
            uiwait(tool.Dialog)
        end
        
        function btnDelUnder(tool,~,~)
            if isempty(tool.DelUnder)
                disp('Type row number to erase.');
                return;
            end
            A=tool.LabelMasks;
            row=str2double(tool.DelUnder.String);
            A(row:end,:)=0;
            tool.LabelMasks=A;
            
        end
        
        function keyPressCB(tool,~,event)
            switch event.Key
                case 's'
                    if tool.RadioErase.Value==1
                        tool.RadioDraw.Value=1;
                        tool.RadioErase.Value=0;
                    elseif tool.RadioDraw.Value==1
                        tool.RadioDraw.Value=0;
                        tool.RadioErase.Value=1;
                    end
                case {'1','2','3','4','5','6','7','8','9'}    
                    tool.PenSize=str2double(event.Key);
                    set(tool.Slider,'Value',str2double(event.Key));
                    set(tool.TransparencyHandle, 'AlphaData', 0.3*tool.LabelMasks);
                case 't'
                    tool.btnTESTPush;
                case 'b'
                    tool.btnBotPush;
                    
            end
        end

        function btnBotPush(tool,~,~)
            Msk=tool.LabelMasks>0.1;
            if isempty(find(Msk, 1))
                return;
            end
            [Msk]=extendBottomForNext20(tool.LabelMasksPrev,Msk,str2double(tool.DelUnder.String));
            tool.LabelMasks=Msk;
            tool.btnTESTPush;
        end
        
        function btnTESTPush(tool,~,~)
            Msk=tool.LabelMasks>0.1;
            if isempty(find(Msk, 1))
                return;
            end
            UMsk=false(size(tool.LabelMasks));
            [Length,Msk,IdxList]=IAB_Msk2Length(Msk,tool.minAL);%meanwhile, small things are removed.
            [lkg,lng,Ulist,UMsk,Msk,PoleSheet]=IAB_Length2LNG(Length,UMsk,IdxList,Msk);%meanwhile, top divLost cases are removed.
            
            [Msk,modFrms]=modChangedCell(tool.LabelMasksPrev,Msk,str2double(tool.DelUnder.String));
            %<tool.LabelMasksPrev> is old mask [double] 0 and 1
            %<Msk> is the new mask [logical]
            %<str2double(tool.DelUnder.String)>is delUnderRow.
            %Aim: Find out which cellID is changed manually and morphologically modify its shape.
            tool.LabelMasks=double(Msk);
            set(tool.U1,'String',num2str(size(Ulist,1))); 
            Ulist2=tool.prevUlist2;
            if size(Ulist,1)==0 %no cell needs manual correction, mark all isolated cells. Create lineages but do not save (leave the saving step to 'Done').
                [PHparts,PHMap]=IAB_simpDecipherLineagev2(lng);
                if tool.isWasInUlist
                    [Ulist2,ICv2,UMsk2]=IAB_testFittingErrv2(Msk,PHparts);%test fitting problems for all frames. Create ICv2 for the first time.
                else
                    if isempty(modFrms) && ~isempty(Ulist2)
                        disp('Pressed T twice.');return;
                    else
                        [Ulist2,ICv2,UMsk2]=IAB_testFittingErrSelFrms(Msk,PHparts,modFrms,tool.prevICv2,tool.prevUlist2,tool.prevUmsk2);%test fitting problems for selected frames
                    end
                end
                tool.prevUlist2=Ulist2;
                tool.prevUmsk2=UMsk2;
                tool.prevICv2=ICv2;
                
                if isempty(Ulist2)%All good. Proceed to done.
                    if size(PHMap,1)==1
                        J=cat(3,double(UMsk*255),zeros(size(UMsk)),ones(size(UMsk)));
                        set(tool.TransparencyHandle, 'CData',J);
                        set(tool.TransparencyHandle, 'AlphaData', 0.3*tool.LabelMasks);
                        return;
                    end
                    [DMsk,coor4Del,ParCoor,isParIso]=IAB_drawGreenMask(size(Msk),PHparts,PHMap,lkg,IdxList);
                    [LineISO,LineLNG]=IAB_createNewLineBatch(size(Msk),ParCoor,isParIso);%create lineage line for div cycles, (:,:,1) is isolated cells, (:,:,2) is lineaged cells.
                    tool.LineageToSave={lkg,lng,single(PoleSheet),coor4Del,ICv2};% pass things to be saved in lkg.mat to tool.LineageToSave
                    J=cat(3,zeros(size(DMsk)),double((DMsk|LineISO)*255),double((Msk|LineLNG)*255));% the green cells are not in complete division cycles.
                    AlphaX=0.3*tool.LabelMasks+0.2*(LineISO | LineLNG);
                    set(tool.TransparencyHandle, 'CData',J);
                    set(tool.TransparencyHandle, 'AlphaData', AlphaX);
                else
                    %if Ulist2 is not empty, create new Umsk:
                    J=cat(3,double(UMsk2*255),zeros(size(UMsk2)),ones(size(UMsk2)));
                    set(tool.TransparencyHandle, 'CData',J);
                    set(tool.TransparencyHandle, 'AlphaData', 0.3*tool.LabelMasks);

                end
                tool.isWasInUlist=false;
            else % Ulist is not empty. Still need manual correction.
                J=cat(3,double(UMsk*255),zeros(size(UMsk)),ones(size(UMsk)));
                set(tool.TransparencyHandle, 'CData',J);
                set(tool.TransparencyHandle, 'AlphaData', 0.3*tool.LabelMasks);
                tool.isWasInUlist=true;
                tool.prevUlist2=[];
                tool.prevICv2=[];
            end

            if ~isempty(Ulist)
                PixA=Ulist(1,1);
                api = iptgetapi(tool.Scroll);
                api.setMagnificationAndCenter(3,PixA*32,128);%<-------------set magnification and center.
            end
            if ~isempty(Ulist2)
                PixA=Ulist2(1,1);
                api = iptgetapi(tool.Scroll);
                api.setMagnificationAndCenter(3,PixA*32,128);%<-------------set magnification and center.
                set(tool.U1,'String',['[2]',num2str(size(Ulist2,1))]); 
            end


            tool.LabelMasksPrev=double(Msk);
        end
        
        function btnEmptizePush(tool,~,~)
            tool.LabelMasks=false(size(tool.LabelMasks));
            tool.U1.String='0';
        end
        
        function btnDelAfter(tool,~,~)
            DANum=str2double(tool.DelAfter.String);
            pale=zeros(size(tool.LabelMasks));
            pale(:,1:DANum*32)=tool.LabelMasks(:,1:DANum*32);
            tool.LabelMasks=pale;
            set(tool.TransparencyHandle, 'AlphaData',0.3*tool.LabelMasks);
        end
        
        function changeSliderRange(tool,src,~)%<----define pen size
            value = str2double(src.String);
            if strcmp(src.Tag,'sliderMin')
                tool.Slider.Min = value;
                tool.Slider.Value = value;
                tool.PenSize = value;
                tool.PenSizeText.String = sprintf('%d',value);
            elseif strcmp(src.Tag,'sliderMax')
                tool.Slider.Max = value;
                tool.Slider.Value = value;
                tool.PenSize = value;
                tool.PenSizeText.String = sprintf('%d',value);
            end
        end
        
        function radioDraw(tool,src,~)%<----keep
            tool.RadioErase.Value = 1-src.Value;
        end
        
        function radioErase(tool,src,~)%<----keep
            tool.RadioDraw.Value = 1-src.Value;
        end
        
        function continuousSliderManage(tool,~,callbackdata)%<----keep
            tag = callbackdata.AffectedObject.Tag;
            if strcmp(tag,'pss')%<--------------if it is the pen size slider
                tool.PenSize = round(callbackdata.AffectedObject.Value);
                ps = tool.PenSize;
                tool.PenSizeText.String = sprintf('%d',ps);
                [Y,X] = meshgrid(-ps:ps,-ps:ps);
                Mask = sqrt(X.^2+Y.^2) < ps;
                r1 = ceil(tool.Axis.YLim(1));
                r2 = floor(tool.Axis.YLim(2));
                c1 = ceil(tool.Axis.XLim(1));
                c2 = floor(tool.Axis.XLim(2));
                rM = round(mean(tool.Axis.YLim));
                cM = round(mean(tool.Axis.XLim));
                tool.TransparencyHandle.AlphaData(max(1,r1):min(size(tool.Image,1),r2),max(1,c1):min(size(tool.Image,2),c2)) = 0;
                if r1 >= 1 && r2 <= size(tool.Image,1) && c1 >= 1 && c2 <= size(tool.Image,2) ...
                        && rM-ps >= 1 && rM+ps <= size(tool.Image,1) && cM-ps >=1 && cM+ps <= size(tool.Image,2)
                    tool.TransparencyHandle.AlphaData(rM-ps:rM+ps,cM-ps:cM+ps) = Mask;
                end
            else%<-----------------------------if it is the threshold slider(disabled already)
                value = callbackdata.AffectedObject.Value;
                if strcmp(tag,'uts')
                    tool.UpperThreshold = value;
                elseif strcmp(tag,'lts')
                    tool.LowerThreshold = value;
                end
                I = tool.Image;
                I(I < tool.LowerThreshold) = tool.LowerThreshold;
                I(I > tool.UpperThreshold) = tool.UpperThreshold;
                I = I-min(I(:));
                I = I/max(I(:));
                tool.ImageHandle.CData = I;
            end
        end
        
        function sliderManage(tool,src,~)%<----keep
            tool.PenSize = round(src.Value);
            set(tool.TransparencyHandle, 'AlphaData', 0.3*tool.LabelMasks);
        end
 
        function popupManage(tool,src,~)%<----keep
            tool.LabelIndex = src.Value;
            set(tool.TransparencyHandle, 'AlphaData', 0.3*tool.LabelMasks);
        end
        
        function closeDialog(tool,~,~)%<----keep
            delete(tool.Dialog);
            delete(tool.Figure);
        end
        
        function buttonDonePushed(tool,~,~)%<-important function.
            tool.btnTESTPush;
            if ~strcmp(tool.U1.String,'0')
                return;
            end
            NoOverlap = sum(tool.LabelMasks,3) <= 1;
            tool.LabelMasks = tool.LabelMasks.*repmat(NoOverlap,[1 1 tool.NLabels]);
            tool.DidAnnotate = 1;
            delete(tool.Dialog);
            delete(tool.Figure);
        end
        
        function btnRstPushed(tool,~,~)
            tool.LabelMasks=tool.RawMsk;
            set(tool.TransparencyHandle, 'AlphaData',0.3*tool.LabelMasks);
            
        end
        
        function closeFigure(tool,~,~)%<----keep
            delete(tool.Figure);
            delete(tool.Dialog);
        end
        
        function mouseMove(tool,~,~)%<-important function.
            p = tool.Axis.CurrentPoint;%<-current coor of pointer.
            col = round(p(1,1));
            row = round(p(1,2));
            set(tool.PX,'String',num2str(col));
            set(tool.PY,'String',num2str(row));
           
            ps = tool.PenSize;
            if row > ps && row <= size(tool.Image,1)-ps && col > ps && col <= size(tool.Image,2)-ps && tool.MouseIsDown
                [Y,X] = meshgrid(-ps:ps,-ps:ps);
                Curr = tool.LabelMasks(row-ps:row+ps,col-ps:col+ps,tool.LabelIndex);
                Mask = sqrt(X.^2+Y.^2) < ps;
                if tool.RadioDraw.Value == 1
                    tool.LabelMasks(row-ps:row+ps,col-ps:col+ps,tool.LabelIndex) = max(Curr,Mask);%<----OR operation
                    tool.TransparencyHandle.AlphaData(row-ps:row+ps,col-ps:col+ps) = 0.3*max(Curr,Mask);
                elseif tool.RadioErase.Value == 1
                    tool.LabelMasks(row-ps:row+ps,col-ps:col+ps,tool.LabelIndex) = min(Curr,1-Mask);
                    tool.TransparencyHandle.AlphaData(row-ps:row+ps,col-ps:col+ps) = min(Curr,0.3*(1-Mask));
                end
            end
        end
        
        function mouseDown(tool,~,~)%<----keep
            tool.MouseIsDown = true;
        end
        
        function mouseUp(tool,~,~)%<----keep
            tool.MouseIsDown = false;
        end
    end
end
