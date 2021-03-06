
function bwpaint1(A,original)

    figure(1), subplot(121), imshow(original, []), colormap(jet);

    subplot(122), imagesc(A);

    colormap('default'); 
    nz = 1;
    if nz>1
        h0 = uicontrol('Style', 'text', ...
             'Position', [100 0 200 50], ...
             'String','Slice No. (slide to change image slice',...
             'Tag','txtErase');
        h = uicontrol('Style', 'slider', ...
             'Position', [100 50 200 50], ...
             'Value',0,...
             'Max',nz-1,...
             'Min',0,...
             'Tag','sliderZ',...
             'SliderStep',[(1/(nz-1)) (1/(nz-1))],...
             'Callback', @display);
    else
        h = uicontrol('Style', 'slider', ...
         'Position', [100 50 200 50], ...
         'Value',0,...
         'Max',0,...
         'Min',0,...
         'Tag','sliderZ',...
         'SliderStep',[1 1],...
         'Callback', @display);
    end
    h1 = uicontrol('Style', 'pushbutton', ...
         'Position', [100 100 200 50], ...
         'String','Export image to A.mat...',...
         'Callback', @export);
    h4 = uicontrol('Style', 'text', ...
         'Position', [100 250 200 50], ...
         'String','Brush Size (slide changes brush size, left click on image to paint, right click to erase...',...
         'Tag','txtErase');

    h5 = uicontrol('Style', 'slider', ...
         'Position', [100 300 200 50], ...
         'Value',0,...
         'Max',5,...
         'Min',0,...
         'Tag','sizeErase',...
         'SliderStep',[0.2  0.2]);


    axis image;
    m=size(A,1);n=size(A,2);%size of the image

    xlim([1 n]);
    ylim([1 m]);


    % Unpack gui object
    gui = get(gcf,'UserData');
    set(gcf,'Pointer','arrow');

    % Make a fresh figure window
    set(gcf,'WindowButtonDownFcn',@startmovit);

    % Store gui object
    set(gcf,'UserData',{gui;A});

function export(src,evnt)
    temp = get(gcf,'UserData');
    gui=temp{1};
    A=temp{2};
    manual_paint = A;
    save('manual_paint.mat','manual_paint');

function display(src,evnt)
    temp = get(gcf,'UserData');
    gui=temp{1};
    A=temp{2};
    a=findobj(gcf,'Style','Slider','-and','Tag','sliderZ');
    % imagesc(A(:,:,round(get(a,'Value'))*3+1:round(get(a,'Value'))*3+3));
    imagesc(A);
    colormap('default');
    axis image;

function startmovit(src,evnt)
    temp = get(gcf,'UserData');
    gui=temp{1};
    A=temp{2};

    hr=findobj(gcf,'Style','Slider','-and','Tag','sizeErase');
    r=get(hr,'Value');
    a=findobj(gcf,'Style','Slider','-and','Tag','sliderZ');
    im2 = A;

    pos = get(gca,'CurrentPoint');
    flag_btn=get(src,'SelectionType');

    [m,n]=size(A);

    cm=round(pos(1,2));
    cn=round(pos(1,1));

    if strcmp(flag_btn,'normal')
    im2(max(cm-r,1):min(cm+r,m),max(cn-r,1):min(cn+r,n))=1; 
    elseif strcmp(flag_btn,'alt')
    im2(max(cm-r,1):min(cm+r,m),max(cn-r,1):min(cn+r,n))=0; 
    end
    
    A = im2;
    imagesc(A);
    colormap('default');
    axis image;

% Set callbacks
    gui.currenthandle = src;
    thisfig = gcbf();
    set(thisfig,'WindowButtonMotionFcn',@movit);
    set(thisfig,'WindowButtonUpFcn',@stopmovit);
 
% Store gui object
    set(gcf,'UserData',{gui;A});



function movit(src,evnt)
% Unpack gui object
    ud = get(gcf,'UserData');
    flag_btn=get(src,'SelectionType');
    try
        if isequal(gui.startpoint,[])
            return
        end
    catch
    end
 
    pos = get(gca,'CurrentPoint');
    img=ud{2};
    hr=findobj(gcf,'Style','Slider','-and','Tag','sizeErase');
    r=get(hr,'Value');
    a=findobj(gcf,'Style','Slider','-and','Tag','sliderZ');
    % im2=img(:,:,round(get(a,'Value'))*3+2);
    im2 = img;
    [m,n]=size(im2);

    cm=round(pos(1,2));
    cn=round(pos(1,1));

    if strcmp(flag_btn,'normal')
    im2(max(cm-r,1):min(cm+r,m),max(cn-r,1):min(cn+r,n))=1; 
    elseif strcmp(flag_btn,'alt')
    im2(max(cm-r,1):min(cm+r,m),max(cn-r,1):min(cn+r,n))=0; 
    end
    
    img = im2;
    imagesc(img);
    colormap('default');
    axis image;

    ud{2}=img;
 
% Store gui object
    set(gcf,'UserData',ud);
 

function stopmovit(src,evnt)
% Clean up the evidence ...
    thisfig = gcbf();
    ud = get(gcf,'UserData'); 
    set(thisfig,'WindowButtonUpFcn','');
    set(thisfig,'WindowButtonMotionFcn','');
    axis image;
    set(gcf,'UserData',ud);