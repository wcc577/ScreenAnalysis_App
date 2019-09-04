function Screen_Analysis(app)
% Screen_CSV.m Screen .csv data caculation of SCR Brightness 2017/10/30 jwu
% Adding 9 points 2D optical info    ver.180410
% Support .mat file format  input 181210


n2f=0.291863508; % conversion 1Nits to FL
nv=[95; 90; 80; 70; 50; 30; 5]; nv=nv/100; % contour level
Sr=0.1; % shrinkage ROI for finding max peak. 0.1,0.1 --> (16mm)0.0
% app.smo=60;%30 % Gaussian smoothing number 5,10-->(16mm)10

% app.NP9=0; % import 9 points optical data
Rm=1;  % Rm=1 setup the ROI for each 9 point view
%readF=1; % read file or not 1=Yes
ssL=1;   % SCR ROI select: 0: No. 1=Right. 2=Left. 3=R+L.
sa=0; % fixed Scale of SCR if NOT 0
dotp=0; % peak SCR print with dot


if isempty(app.file2D)==0
    fmt=app.file2D(end-2:end); 
    file=[app.fileDir app.file2D];
    if fmt=='csv'
        B_2D=round(csvread(file,2,1)*n2f,1);
    elseif fmt=='mat'
        load(file);
        B_2D=round(data*n2f,1); 
    end

    [m n]=size(B_2D); ROI=[1; m; 1; n]; % Default full range [Y-min; Y-max; X-min; X-max]

    figure(1)
    imagesc(B_2D); colormap 'jet'
    title('Click the Up-Left and Bottom-Right points for 2D ROI select');
    [x,y]=ginput(2); x=round(x); y=round(y); ROI=[min(y); max(y); min(x); max(x)];

    B_2D=B_2D(ROI(1):ROI(2),ROI(3):ROI(4)); 
    B_2D=smoothdata(B_2D,'gaussian',app.smo);

    [ms ns]=size(B_2D);
    ROIs=round([ms*(Sr); ms*(1-Sr); ns*(Sr); ns*(1-Sr)]); ROIs(find(ROIs==0))=1;
    ROIsM=zeros(ms, ns); ROIsM(ROIs(1):ROIs(2),ROIs(3):ROIs(4))=1;
end

%  if noF>2
%     if readF==1,
if isempty(app.file3DRB)==0
    fmt=app.file3DRB(end-2:end); 
    file=[app.fileDir app.file3DRB];
    if fmt=='csv'
        R_Bright=round(csvread(file,2,1)*n2f,1);
    elseif fmt=='mat'
        load(file);
        R_Bright=round(data*n2f,1); 
    end
end

if isempty(app.file3DRD)==0
    fmt=app.file3DRD(end-2:end); 
    file=[app.fileDir app.file3DRD];
    if fmt=='csv'
        R_Dark=round(csvread(file,2,1)*n2f,1);
    elseif fmt=='mat'
        load(file);
        R_Dark=round(data*n2f,1); 
    end
end

if isempty(app.file3DLB)==0
    fmt=app.file3DLB(end-2:end); 
    file=[app.fileDir app.file3DLB];
    if fmt=='csv'
        L_Bright=round(csvread(file,2,1)*n2f,1);
    elseif fmt=='mat'
        load(file);
        L_Bright=round(data*n2f,1); 
    end
end

if isempty(app.file3DLD)==0
    fmt=app.file3DLD(end-2:end); 
    file=[app.fileDir app.file3DLD];
    if fmt=='csv'
        L_Dark=round(csvread(file,2,1)*n2f,1);
    elseif fmt=='mat'
        load(file);
        L_Dark=round(data*n2f,1); 
    end
end


    if ssL==1||3        
        figure(9)
        imagesc(R_Bright); colormap 'jet'
        title('Click the Up-Left and Bottom-Right pos for Right-Eye ROI select');
        [x,y]=ginput(2); x=round(x); y=round(y); ROI=[min(y); max(y); min(x); max(x)];
        R_Bright=R_Bright(ROI(1):ROI(2),ROI(3):ROI(4));
        R_Dark=R_Dark(ROI(1):ROI(2),ROI(3):ROI(4));

        [ms, ns]=size(R_Bright);
        ROIs=round([ms*(Sr); ms*(1-Sr); ns*(Sr); ns*(1-Sr)]); ROIs(find(ROIs==0))=1;
        ROIsMR=zeros(ms, ns); ROIsMR(ROIs(1):ROIs(2),ROIs(3):ROIs(4))=1;
    else
        ROIsMR=ROIsM; 
        R_Bright=R_Bright(ROI(1):ROI(2),ROI(3):ROI(4));       
        R_Dark=R_Dark(ROI(1):ROI(2),ROI(3):ROI(4));
    end
    
    if ssL==2||3
        figure(9)
        imagesc(L_Bright); colormap 'jet'
        title('Click the Up-Left and Bottom-Right pos for Left-Eye ROI select');
        [x,y]=ginput(2); x=round(x); y=round(y); ROI=[min(y); max(y); min(x); max(x)];
        L_Bright=L_Bright(ROI(1):ROI(2),ROI(3):ROI(4));
        L_Dark=L_Dark(ROI(1):ROI(2),ROI(3):ROI(4));

        [ms ns]=size(L_Bright);
        ROIs=round([ms*(Sr); ms*(1-Sr); ns*(Sr); ns*(1-Sr)]); ROIs(find(ROIs==0))=1;
        ROIsML=zeros(ms, ns); ROIsML(ROIs(1):ROIs(2),ROIs(3):ROIs(4))=1;
    else
        ROIsML=ROIsM;
        L_Bright=L_Bright(ROI(1):ROI(2),ROI(3):ROI(4));
        L_Dark=L_Dark(ROI(1):ROI(2),ROI(3):ROI(4));
    end
    
    R_Bright=smoothdata(R_Bright,'gaussian',app.smo); %adding smooth 190824
    L_Bright=smoothdata(L_Bright,'gaussian',app.smo);
    R_Darkt=smoothdata(R_Dark,'gaussian',app.smo); %adding smooth 190824
    L_Dark=smoothdata(L_Dark,'gaussian',app.smo);

    L_SCR=round(L_Bright./L_Dark); 
    A=find(L_Dark==0); L_SCR(A)=0;  %remove Inf
    L_SCR=smoothdata(L_SCR,'gaussian',app.smo);

    R_SCR=round(R_Bright./R_Dark); 
    B=find(R_Dark==0); R_SCR(B)=0;  %remove Inf
    R_SCR=smoothdata(R_SCR,'gaussian',app.smo);
% end

fg1=figure(1); clf; fg1.Position=[560 528 560 420];
colormap 'parula'
% h=fspecial('gaussian',sigma,sigma);
% dataFiltered=imfilter(B_2D,h,'replicate','same','conv');
imagesc(B_2D); hold on;  %colormap 'jet';
s2DB=max(max(B_2D)); s=round(s2DB,1); nc=round(nv*s);
[C,h]=contour(B_2D,'ShowText','on','LineColor',[0.9 0.9 0.9],'LabelSpacing',800); clabel(C,h); 
til=['2D Brightness of ' app.Auditorium ' (Max. Brightness '...
    num2str(s) ' FL)'];
title(til); hold off; axis('equal'); xlabel('X Axis');
%fg1.Position=[560 528 560 420];
saveas(gcf,[app.fileDir '2D Brightness.png'])

noF=4;
if noF>2
    fg2=figure(2); fg2.Position=[1120 120 560 820];
    subplot(3,1,1)
    imagesc(R_Bright); hold on
    sRB=max(max(R_Bright.*ROIsMR)); s=round(sRB,1);
    [C,h]=contour(R_Bright,'ShowText','on','LineColor',[0.9 0.9 0.9],'LabelSpacing',800); %clabel(C,h); 
    til=['Right-Eye 3D Brightness of ' app.Auditorium ' (Max. Brightness ' ...
        num2str(s) ' FL)'];
    title(til); hold off; axis('equal'); caxis([0 s]);

    subplot(3,1,2)
    imagesc(R_SCR); hold on 
    sRS=max(max(R_SCR.*ROIsMR)); s=round(sRS); 
    if sa~0; s=sa; end %%%% nc level on Right-eye
    nc=round(nv*s);  
    [C,h]=contour(R_SCR,nc,'ShowText','on','LineColor',[0.9 0.9 0.9],'LabelSpacing',800); %clabel(C,h); 
    til2=['Right-Eye 3D SCR of ' app.Auditorium ' (Max. SCR ' num2str(s) ':1)'];
    title(til2); hold off; axis('equal'); caxis([0 s]);

    subplot(3,1,3)%, title('Cross Peak SCR')
    [a b]=max(R_SCR.*ROIsMR,[],2); %a2=smooth(a,0.05,'rloess'); 
    [sRS c]=max(a); 
    if dotp==1;
        yyaxis right; plot(b,a,'ro'),ylabel('Peak SCR')
        yyaxis left; 
    end
    plot(R_SCR(c,:)),ylabel('Cross SCR'),grid on
    xlim([0 ns]); a=ylim; ylim([0 a(2)]);
    til3=['Right-Eye Peak 3D SCR Crossection of ' app.Auditorium];
    title(til3);
    %fg2.Position=[1 1 560 820];
    saveas(gcf,[app.fileDir '-3D Right-Eye SCR.png'])

    fg3=figure(3); fg3.Position=[25 120 560 820];
    subplot(3,1,1)
    imagesc(L_Bright); hold on; %colormap 'jet';
    sLB=max(max(L_Bright.*ROIsML)); s=round(sLB,1);
    [C,h]=contour(L_Bright,'ShowText','on','LineColor',[0.9 0.9 0.9],'LabelSpacing',800); %clabel(C,h); 
    til=['Left-Eye 3D Brightness of ' app.Auditorium ' (Max. Brightness ' ...
        num2str(s) ' FL)'];
    title(til); hold off; axis('equal'); caxis([0 s]);

    subplot(3,1,2)
    imagesc(L_SCR); hold on
    sLS=max(max(L_SCR.*ROIsML)); s=round(sLS); 
    if sa~0; s=sa; end %%%% nc level on Right-eye
    nc=round(nv*s);
    [C,h]=contour(L_SCR,nc,'ShowText','on','LineColor',[0.9 0.9 0.9],'LabelSpacing',800); %clabel(C,h);
    til2=['Left-Eye 3D SCR of ' app.Auditorium ' (Max. SCR ' num2str(s) ':1)'];
    title(til2); hold off; axis('equal'); caxis([0 s]);

    subplot(3,1,3)%, title('Cross Peak SCR')
    [a b]=max(L_SCR.*ROIsML,[],2); %a2=smooth(a,0.05,'rloess'); 
    [sLS c]=max(a); 
    if dotp==1;
        yyaxis right; plot(b,a,'ro'),ylabel('Peak SCR')
        yyaxis left; 
    end
    plot(L_SCR(c,:)),ylabel('Cross SCR'),grid on
    xlim([0 ns]); a=ylim; ylim([0 a(2)]);
    til3=['Left-Eye Peak 3D SCR Crossection of ' app.Auditorium];
    title(til3);
    %fg3.Position=[25 120 560 820];
    saveas(gcf,[app.fileDir '-3D Left-Eye SCR.png'])
end
app.MSG{1}=sprintf('   Theater:  %s \n',app.Auditorium);
app.MSG{2}=sprintf('   Seat: with Westboro P230U_16 8mm lens \n'];
app.MSG{3}=sprintf('   2D Peak Brightness: %2.1f FL\n',s2DB];
if noF>2
    app.MSG{4}=sprintf('   3D Peak Brightness Left/Right Eye: %2.1f/%2.1f FL \n',sLB,sRB);
    app.MSG{5}=sprintf('   3D Peak SCR Left/Right Eye: %3.0f / %3.0f :1 \n',sLS,sRS);
end

    fid=fopen([app.fileDir 'Summary.txt'],'w');
    fprintf(fid,' ...File name %s ...\n\n',file);
    fprintf(fid,'   Theater: %s \n',app.Auditorium);
    fprintf(fid,'   Seat: with Westboro P230U_16 8mm lens \n');
    fprintf(fid,'   2D Peak Brightness: %2.1f FL\n',s2DB);
    if noF>2
        fprintf(fid,'   3D Peak Brightness Left/Right Eye: %2.1f/%2.1f FL \n',sLB,sRB);
        fprintf(fid,'   3D Peak SCR Left/Right Eye: %3.0f / %3.0f :1 \n',sLS,sRS);
    end
    fclose(fid);

    if app.NP9==1;
        n9=length(P9N);
        %if app.NP9==1
            fg9=figure(9); clf;
            fg9.Position=[560 120 1025 820];
            for j=1:n9
                file=[app.fileDir P9N{j}];

                if readF==1;            
                    %ps9=smoothdata(ps9,'gaussian',smo);           
                    if fmt=='csv'
                        ps9=csvread(file,2,1)*n2f; 
                    elseif fmt=='mat'
                        load(file);
                        ps9=data*n2f; 
                    end            
                else
                    ps9=data(:,:,Id(j+5))*n2f;
                end

                if Rm==1 
                    imagesc(ps9);
                    title('Click the Up-Left and Bottom-Righ points for ROI select');
                    [x,y]=ginput(2); x=round(x); y=round(y); 
                    [m n]=size(ps9);
                    A=[min(y); max(y); min(x); max(x)];
                    if A(1)<1; A(1)=1; end
                    if A(3)<1; A(3)=1; end
                    if A(2)>m; A(2)=m; end
                    if A(4)>n; A(4)=n; end
                    MROI(:,j)=A;           
                else
                    MROI(:,j)=[1; m; 1; n];
                end
                P9Bright{j}=ps9(MROI(1,j):MROI(2,j),MROI(3,j):MROI(4,j));

            end

            for j=1:n9
                subplot(3,3,j)
                imagesc(P9Bright{j}); hold on; 
                P9PB(j)=max(max(P9Bright{j})); 
                [C,h]=contour(P9Bright{j},'ShowText','on','LineColor',[0.9 0.9 0.9],'LabelSpacing',300); clabel(C,h);
                hold off
                til=P9N{j}; 
                if readF==1; til=til(1:end-4); end
                til=[': ' til];
                if j==2; til=['9 Points 2D Brightness of ' app.Auditorium '/' til];  end
                title(til);
                %axis('equal')
            end
            %fg9.Position=[560 120 1025 820];
            saveas(gcf,[app.fileDir '9 Points 2D Brightness.png'])
        %end
    end
end