function rotAngle=rotationGUI(I)
    %# read image
%     I = imread('cameraman.tif');

    %# setup GUI
    hFig = figure('menu','none'); title('Rotate Image Angle');
    hAx = axes('Parent',hFig);
    hSld=uicontrol('Parent',hFig, 'Style','slider', 'Value',0, 'Min',-5,...
        'Max',5, 'SliderStep',[0.02 0.02], ...
        'Position',[150 5 300 20], 'Callback',@slider_callback); 
    hTxt = uicontrol('Style','text', 'Position',[290 28 20 15], 'String','0');
    hBtn = uicontrol('Style', 'pushbutton', 'String', 'Done',...
        'Position', [20 20 50 20],...
        'Callback', @button_callback);

    %# show image
    rotAngle=0;
    imagesc(I, 'Parent',hAx)
    

    waitfor(hFig)
    
    %# Callback function
    function slider_callback(hObj, eventdata)
        rotAngle = get(hObj,'Value');        %# get rotation angle in degrees
        imagesc(imrotate(I,rotAngle), 'Parent',hAx)  %# rotate image
        set(hTxt, 'String',num2str(rotAngle))       %# update text
%         assignin('base', 'rotAngle', rotAngle)
        
    end

    function button_callback(hObj,eventdata)
        close(gcf);
        return
    end
end