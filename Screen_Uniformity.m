function Screen_Uniformity(app)
    % Screen Uniformity Analysis Screen_Analy.m
    % Read the westboro .mat or .csv file and choice the location for analysis

    % 2017/11/ jwu
    
    n2f=0.291863508; % conversion 1Nits to FL
    avgn=1; %"1" use range averaged for H profile, "0" for NoN (raw data)
    scn=-10; %scan pixel pitch, '0' for NoN, <0 means average; ONLY works when avgn=0;
    sigma=60;%120(R16)&for Home cinema "it is 60 before" 80(VHG); % sigma of gaussian filter for projector profile remove 60, 30 is too small, but off-axis is 30
    gau=0; % Data gaussian filtering, min. is 1, '0' for NoN
    roi=1; % Select ROI=1, 0=input data.
    sfig=1; % save figure or not, 1=Yes

    fmt=0; % fmt==1, usr uniformity for brightness
    readF=1; % read file or not 1=Yes

    imrot=2; % 0= No, 1=fixed angle, 2=select

    if readF==1
        Ro=exist('fileDir');
        if Ro==1;
            [fileName,fileDir]=uigetfile('*.mat;*.csv','Choose the file you want to process.',fileDir);
        else
            [fileName,fileDir]=uigetfile('*.mat;*.csv','Choose the file you want to process.');
        end
        ftyp=fileName(end-2:end);
        switch ftyp
            case 'mat'
                load([fileDir fileName]); datas=data;
            case 'csv'
                 datas=csvread([fileDir fileName],2,1);
        end
    else
        datas=data(:,:,1);
    end
    datas=datas*n2f; %data=data.';

    if imrot==1; %use secified rotating angle
        rotAngle=17.82;%16.4;%17.82;
        datas=imrotate(datas,rotAngle);
    elseif imrot==2; %define rotating angle
        rotAngle=rotationGUI(datas);
        datas=imrotate(datas,rotAngle);
        fprintf('Image Rotating Angle is %3.2f degrees \n',rotAngle)
    else
        rotAngle=0;
    end

    %% ROI select
    if roi==1  %%%%%%%%%%<<< roi control   >>>
        figure(4)
        %imagesc(data,autoScale(0.02, 0.98, data)); colormap gray;
        imagesc(datas);    
        title('Click the Up-Left and Bottom-Righ points for ROI secect');
        [x,y]=ginput(2); x=round(x); y=round(y);
    else
        %x=[544;1769]; y=[50;993]; %R97
        %x=[433;1765]; y=[54;996]; %R16
        %x=[362;1765]; y=[57;996]; %R150X
        %x=[270;1500]; y=[479;1004]; %R160
        %x=[266;1500]; y=[479;1004]; %R160
        %x=[270;1500]; y=[479;1004]; %On axis
        %x=[560;1765]; y=[55;995]; %Off axis
        x=[195;1762]; y=[319;1202]; %Enlarge ON-axis
        %x=[504;1773]; y=[12;1024]; %Enlarge OFF-axis
    end

    [m,n]=size(datas);
    A=[min(y); max(y); min(x); max(x)];
    if A(1)<1; A(1)=1; end
    if A(3)<1; A(3)=1; end
    if A(2)>m; A(2)=m; end
    if A(4)>n; A(4)=n; end
    sm=input('Input the seam number[Enter for 0]: ');
    if isempty(sm)==1; sm=0; end
    %ROI=[1;1;1927;1148];
    %ROI=[375;700;1560;1100];  % ROI of [Y1,X1, Y2,X2]
    datas=datas(A(1):A(2),A(3):A(4));

    [m,n]=size(datas);
    rawData=datas;

    if gau>=1
        datas=imgaussfilt(rawData, gau);
    end

    h=fspecial('gaussian',sigma,sigma);
    %dataFiltered2=imgaussfilt(data,sigma); %edge is not enhanced ng
    dataFiltered=imfilter(datas,h,'replicate','same','conv');

    if fmt==1;
        dataFiltered=dataFiltered/max(max(dataFiltered))*100;
    end

    figure(4)
    %imagesc(dataFiltered,autoScale(0.02, 0.98, dataFiltered)); colormap 'jet'
    imagesc(dataFiltered); colormap 'jet', hold on
    [C,h]=contour(dataFiltered,'ShowText','on','LineColor',[0.9 0.9 0.9],'LabelSpacing',300); clabel(C,h); hold off
    axis('equal'); title('System Brightness Map')
    if sfig==1; saveas(gcf,[fileDir fileName(1:end-4) ' System Brightness Map.png']); end
    %% Uniformity caculation
    data2=datas./dataFiltered;
    data2=data2/mean2(data2);
    uniformityGauss=(data2-1)*100;
    %% Raw data plot and position selection
    figure(5)
    subplot(2,1,1)
     imagesc(rawData,autoScale(0.02, 0.98, rawData))
     [a b]=caxis; caxis([0,b]);
     colormap gray; %colorbar; 

     if avgn==0;
         title('Click on the position for analysis');
         til2=['Raw data analysis']; 
         if scn>0
             [a center]=ginput(1); center=round(center);
         else
             [a center]=max(max(dataFiltered.'));
             if (m-center)<abs(scn); center=m-abs(scn); end
         end
         yl=[center,center];
     elseif avgn==1
         if roi==1 %%%%%%%%%%%%%<<<<   roi control   >>>
             title('Click on the \color{blue}Up \color{black}and \color{red}Low boundary \color{black}position for analysis');
             [a B]=ginput(2); B=round(B); center=round(mean(B));
         else
             B=[150;350]; center=round(mean(B)); %On axis    <-------------------- SET ROI
             %B=[450;650]; center=650;%round(mean(B)); %Off axis
         end
         til2=['Boundary Averaged data analysis'];
         yl=B*ones(1,2); yl=yl.';
     end
         xl=[1;n]; line(xl,yl);

     subplot(2,1,2)
     if scn>1
        A=rawData((center-scn):scn:(center+scn),:);
        yl2=[0; max(A(2,:))];%yl2=[min(A(2,:)); max(A(2,:))];
     elseif scn<0
         scn2=abs(scn); %cente=center; %440;
        A=mean(rawData((center-scn2):scn2:(center+scn2),:),1);
        yl2=[0; max(A)];%yl2=[min(A); max(A)];
     elseif scn==0
     %else
        A=rawData(center,:); 
        %A=mean(rawData,1);
        %A=mean(rawData(B(1):1:B(2),:),1);  %averaged data
        yl2=[0; max(A)];%yl2=[min(A); max(A)];
     end
     plot(A.'); grid, title(til2);
     xlim([1 n]); ylim(yl2*1.1);
     %R160=A;

     subplot(2,1,1)
     if sm>0
        s=['\color{red}Click on the ' num2str(sm) ' Seam position for drawing'];
        title(s); %disp(s); 
        [px py]=ginput(sm);
        px=round(px); py=round(py); 
     end

     til=[fileName(1:end-4) ' Raw Data '];
     title(til); 

     if sm>0
         subplot(2,1,2)
    % if scn>1
    %     A=rawData((center-scn):scn:(center+scn),:);
    %     yl2=[min(A(2,:)) max(A(2,:))];
    % else
    %     %A=rawData(center,:);
    %     %A=mean(rawData,1);
    %     A=mean(rawData(B(1):1:B(2),:),1);  %averaged data
    %     yl2=[min(A) max(A)];
    % end
     %plot(rawData(center,:)); grid
     %plot(A.'); grid, title(til2);
     %xlim([1 n]);
       xl2=[px px]; xl2=xl2.';
       %%%%line(xl2,yl2,'Color','black','LineStyle',':');

       for j=1:sm; tstr{j}=['Seam' num2str(j)]; end
       ty=ones(sm,1)*round(yl2(2)-(yl(2)-yl(1))*0.9);
       %%%%text(px,ty,tstr);
     end

     if scn>1; legend('-50 Position','Corss line','+50 Position'); end
     if sfig==1; saveas(gcf,[fileDir fileName(1:end-4) ' Uniformity Raw Data.png']); end

    %% Analized Figure 
    figure(6)
    subplot(2,1,1)
     imagesc(uniformityGauss,autoScale(0.02, 0.98, uniformityGauss))

    if sm>0
        s=['\color{red}Click on the ' num2str(sm) ' Seam position for drawing'];
        title(s); %disp(s); 
        [px py]=ginput(sm);
        px=round(px); py=round(py); 
        xl2=[px px]; xl2=xl2.';
    end


     til=[fileName(1:end-4) ' Enhanced Gaussian filtered Uniformity, Sigma ' num2str(sigma) ' /blur ' num2str(gau) ' /rotAngle ' num2str(rotAngle)];
     title(til); 
     colormap gray;  line(xl,yl); %colorbar;
     %[px py]=ginput(sm); px=round(px); py=round(py); yl2=[-10 10];

    subplot(2,1,2)
     if scn>1
        A=uniformityGauss((center-scn):scn:(center+scn),:);
     elseif scn<0
         scn2=abs(scn);
        A=mean(uniformityGauss((center-scn2):1:(center+scn2),:),1);
     elseif scn==0
     %else   
        A=uniformityGauss(center,:);
        %A=mean(uniformityGauss,1);
        %%A=mean(uniformityGauss(B(1):1:B(2),:),1);  %averaged data
     end
     %A=uniformityGauss((center-scn):scn:(center+scn),:);
     plot(A.'-1.5); grid, 
     xlim([1 n]); ylim([-5 5]); %ylim([0 15]);
     yticks([-5 -2.5 -1 0 1 2.5 5]); yticklabels({'-5','-2.5','-1 Low','0','+1 Up','2.5','5'});
     %yticks([0 7.5 10 12.5 15]); yticklabels({'0','7.5','10 Base','12.5','15'});
     yl2=[-5; 5]; title(til2);

     if sm>0
       line(xl2,yl2,'Color','black','LineStyle',':');
       ty=ones(sm,1)*round(max(yl2)*0.9);
       text(px,ty,tstr);
     end

     if scn>1; legend('-50 Position','Corss line','+50 Position'); end
     if sfig==1; saveas(gcf,[fileDir fileName(1:end-4) ' Enhanced Gaussian filtered Uniformity.png']); end

     fprintf('\n Completed with %s \n %s \n',fileDir,fileName)
 
end