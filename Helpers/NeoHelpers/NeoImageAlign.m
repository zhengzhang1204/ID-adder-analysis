classdef NeoImageAlign < handle    
    properties
        Figure
        Dialog
        Axis
        TransparencyHandle
        ImageHandle
        MouseIsDown
        MouseIsUp
        Scroll
        sz
        driftX
        Rn
        rawRng
        newRng
        CrossHair
        PX
        PY
        PXFL
        PYFL
%table data.
        tab
        tabData
        tabSel
    end
    
    methods
        function tool = NeoImageAlign(M,R,pad,rawRng)
            tool.MouseIsDown = false;
            tool.MouseIsUp = true;
            tool.driftX=[0,0];
            tool.sz=size(R);
            tool.rawRng=rawRng;
            tool.newRng=[];
            tool.tabData=[0,0,0,0];
            tool.CrossHair=pad;
            tool.Figure = figure(...%'MenuBar','none', ...
                                 'NumberTitle','off', ...
                                 'Name','Move and selete aligned points', ...
                                 'Position',[330 100 1500 900],...% 'CloseRequestFcn',@tool.closeFigure, ...
                                 'WindowButtonMotionFcn', @tool.mouseMove, ...%click down and hold
                                 'WindowButtonDownFcn', @tool.mouseDown, ...%click down
                                 'WindowKeyPressFcn',@tool.keyPressCB,...
                                 'Resize','on');

            tool.Axis = axes('Parent',tool.Figure,'Position',[0 0 1 1]);
            %How to generate image for display?
            %out=cat(3,RJ+canvas(:,:,1)+Contour,RJ+2*canvas(:,:,2)+canvas(:,:,1),RJ);%RJ:FLTS; canvas:x-cross and text labels; contour:cellPHMaskContour
            %J=cat(3,RJ+Contour,RJ+2*canvas(:,:,2),RJ);
            %So, out=cat(3,tool.Image(:,:,1)+tool.LineLayer,tool.Image(:,:,2)+tool.LineLayer,tool.Image(:,:,3));
            tool.ImageHandle = imshow(imadjust(M),'Parent',tool.Axis);%This stores the image (tool.ImageHandle)
            tool.Scroll = imscrollpanel(tool.Figure,tool.ImageHandle);

            api = iptgetapi(tool.Scroll);
            api.setMagnificationAndCenter(1.5,256,100);%<-------------set magnification and center.
            Ra=imadjust(R);%auto adjust contrast.
            Ra=cat(3,Ra,Ra,zeros(size(Ra)));%<-change color here!
            hold on
            tool.TransparencyHandle = imshow(Ra,'Parent',tool.Axis);%This stores the blue layer [0,0,1] (tool.TransparencyHandle)
            hold off
            set(tool.TransparencyHandle, 'AlphaData', 0.4*ones(tool.sz));
            tool.Rn=Ra;

            tool.Dialog =dialog('WindowStyle', 'normal',...
                                'Name', 'ImFL-Master',...
                                'CloseRequestFcn', @tool.closeDialog,...
                                'Position',[10 50 300 300]);
            tool.tab = uitable('Parent',tool.Dialog,'Data',tool.tabData,'ColumnName',{'Mx','My','FLx','FLy'},'Position',[10 60 220 240],'CellSelectionCallback',@tool.tabSelCB); 
            tool.tab.ColumnWidth={45,45,45,45};
            set(tool.tab,'ColumnEditable',true);

    % XY mouse pointer coordinates.
            uipanel(tool.Dialog,'units','pixels','Position',[235,100,60,80]);
    % PH position:
            uicontrol(tool.Dialog,'Position',[236,168,29,10],'Style','text','HorizontalAlignment','right','String','R.PH:','FontSize',7);
            tool.PY=uicontrol(tool.Dialog,'Position',[264,168, 28,10],'Style','text','HorizontalAlignment','left','String','','FontSize',7);
            uicontrol(tool.Dialog,'Position',[236,156,29,10],'Style','text','HorizontalAlignment','right','String','C.PH:','FontSize',7);
            tool.PX=uicontrol(tool.Dialog,'Position',[264,156, 28,10],'Style','text','HorizontalAlignment','left','String','','FontSize',7);
    % FL position:
            uicontrol(tool.Dialog,'Position',[236,140,29,10],'Style','text','HorizontalAlignment','right','String','R.FL:','FontSize',7);
            tool.PYFL=uicontrol(tool.Dialog,'Position',[264,140, 28,10],'Style','text','HorizontalAlignment','left','String','','FontSize',7);
            uicontrol(tool.Dialog,'Position',[236,128,29,10],'Style','text','HorizontalAlignment','right','String','C.FL:','FontSize',7);
            tool.PXFL=uicontrol(tool.Dialog,'Position',[264,128, 28,10],'Style','text','HorizontalAlignment','left','String','','FontSize',7);

    % pencil/eraser slider
            uicontrol('Parent',tool.Dialog,'Style','text','String','LUT','Position',[10 10 100 20],'HorizontalAlignment','left');
            tool.PenSizeText = uicontrol('Parent',tool.Dialog,'Style','text','String','5','Position',[125 90 50 20],'HorizontalAlignment','center');
            tool.PenSizeText = uicontrol('Parent',tool.Dialog,'Style','text','String','5','Position',[125 90 50 20],'HorizontalAlignment','center');
            uicontrol('Parent',tool.Dialog,'Style','edit','String','10','Position',[240 90 50 20],'HorizontalAlignment','right','Callback',@tool.changeSliderRange,'Tag','sliderMax');
            uicontrol('Parent',tool.Dialog,'Style','edit','String','1','Position',[10 90 50 20],'HorizontalAlignment','left','Callback',@tool.changeSliderRange,'Tag','sliderMin');
            tool.SliderLow = uicontrol('Parent',tool.Dialog,'Style','slider','Min',1,'Max',10,'Value',tool.PenSize,'Position',[10 30 160 20],'Callback',@tool.sliderManage,'Tag','pss');
            tool.SliderHigh = uicontrol('Parent',tool.Dialog,'Style','slider','Min',1,'Max',10,'Value',tool.PenSize,'Position',[10 10 160 20],'Callback',@tool.sliderManage,'Tag','pss');
            addlistener(tool.Slider,'Value','PostSet',@tool.continuousSliderManage);

    % Button Done:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','Done','Position',[250 10 40 40],'Callback',@tool.btnDone);
            % Button Delete spot:
            uicontrol('Parent',tool.Dialog,'Style','pushbutton','String','DelRow','Position',[250 60 40 40],'Callback',@tool.btnDelRow);

            uiwait(tool.Figure)
        end
%% GUI functions:
% Press buttons: AddSpot, Conv2DualSpots, DelSpot, GenTree, GenMap, RefreshPic, Skip, SaveList, SavePic, SLsharp, Test(with TableDone),
        function sliderMove(tool,~,event)

        end




        function keyPressCB(tool,~,event)% place to go if want to set key shortcuts.
        S=tool.Rn;
        drift=tool.driftX;
        ST=get(tool.TransparencyHandle, 'AlphaData');
        switch event.Key
            case 'a' %'a' move left
                S=[S,S(:,1,:)];
                S(:,1,:)=[];
                ST=[ST,ST(:,1)];
                ST(:,1)=[];
                drift=drift+[-1,0];
            case 'd'%'d' move right
                S=[S(:,end,:),S];
                S(:,end,:)=[];
                ST=[ST(:,end),ST];
                ST(:,end)=[];
                drift=drift+[1,0];
            case 'w'%'w' move up
                S=[S;S(1,:,:)];
                S(1,:,:)=[];
                ST=[ST;ST(1,:)];
                ST(1,:)=[];
                drift=drift+[0,-1];
            case 's'%'s' move down
                S=[S(end,:,:);S];
                S(end,:,:)=[];
                ST=[ST(end,:);ST];
                ST(end,:)=[];
                drift=drift+[0,1];
            case 'q'%'q' quit
                
        end
        tool.Rn=S;
        tool.driftX=drift;
        set(tool.TransparencyHandle,'CData',tool.Rn);
        set(tool.TransparencyHandle,'AlphaData',ST);
        end  %done

        function tabSelCB(tool,~,event)
            tool.tabSel=event.Indices;
        end
        function mouseDown(tool,~,~)
            tool.MouseIsDown = true;
            if tool.MouseIsDown && tool.MouseIsUp
                p = tool.Axis.CurrentPoint;
                col = round(p(1,1));
                row = round(p(1,2));
                if row<=0 || col<=0 || row>=tool.sz(1) || col>=tool.sz(2)
                    return;
                end
                uit=[col,row,col-tool.driftX(1),row-tool.driftX(2)];
                if isequal(tool.tabData,[0,0,0,0])
                    tool.tabData=uit;
                else
                    tool.tabData=[tool.tabData;uit];
                end
                set(tool.tab,'Data',tool.tabData);

            % mark crosshair
                if row>4 && col>4 && row<tool.sz(1)-4 && col<tool.sz(2)-4
%                     fullMask=false(tool.sz(1:2));
                    K=get(tool.TransparencyHandle, 'AlphaData');
                    K((row-4):(row+4),(col-4):(col+4))=double(tool.CrossHair);
                    set(tool.TransparencyHandle, 'AlphaData', K);
                end
            end
            tool.MouseIsDown = false;
        end% done

        function mouseUp(tool,~,~)
            tool.MouseIsUp = true;
            tool.MouseIsDown = false;
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

            set(tool.PXFL,'String',num2str(col-tool.driftX(2)));
            set(tool.PYFL,'String',num2str(row-tool.driftX(1)));

        end

        function btnDelRow(tool,src,callbackdata)
            delRow=tool.tabSel(1);
            tool.tabData(delRow,:)=[];
            set(tool.tab,'Data',tool.tabData);
        end
        
        function closeDialog(tool,src,callbackdata)%<----keep
            delete(tool.Dialog);
            delete(tool.Figure);
        end
        
        function btnDone(tool,src,callbackdata)%<----working on
            delete(tool.Dialog);%<-----perserve permanently
            delete(tool.Figure);%<-----perserve permanently
        end
        
        function closeFigure(tool,src,callbackdata)%<----keep
            delete(tool.Figure);
            delete(tool.Dialog);
        end
%% nested functions:


        

    end
end
