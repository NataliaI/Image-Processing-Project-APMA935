function[T, CPUtime]=BSCBInpaint(I, Mask, NSwitches,GetNorm, FigNum)
%--------------------------------------------------------------------------
% BSCBInpaint performs Bertalmio-Sapiro-Caselles-Ballester inpainting
% Alternates between diffusion and image inpainting process
% Calls functions BSCBIter and TVIter for inpainting and diffusion
% Inputs: 
% I: matrix of image to be inpainted 
% Mask: binary mask of inpaining region
% NSwitches: number of times to repeat inpainting and diffusion processes
% GetNorm: =1 if you want to plot the norm of I vs. iterations
% FigNum: figure number for output (string format)
% Outputs: 
% T: max time T (using default step size 0.1)
% CPUtime: CPU time
%--------------------------------------------------------------------------
global dt

% plotting initial image with mask
%----------------------------------
figure;
imagesc(I); axis off; colormap(gray);
title('Initial Image', 'FontSize',12);
saveas(gcf, strcat('BSCBOutput/',FigNum,'InitialImage.eps'),'epsc'); 

MaskVals=I(Mask==1);
Mask_val=MaskVals(1,1);
t=0;

% Applying diffusion to entire image first
%-------------------------------------------
figure; colormap(gray);
disp('Applying 10 iterations of TV to entire image');
lambda_val=100;
I=TVIter(I, zeros(size(Mask)), 10,0,0,t,lambda_val);

% Starting Inpainting
%---------------------
tStart=tic;
t=0; 
A=15; % num iterations of inpainting
B=2;  % num iterations of TV 
if GetNorm==1
    normsI=zeros(1,NSwitches*(A+B)); %keeps track of ||I|| at each iteration
end

for switches=1:NSwitches
    disp(['switch ' num2str(switches) ' of ' num2str(NSwitches)]);
   
    % BSCB Inpainting
    %-----------------
    disp([num2str(A) ' iterations of BSCB inpainting']);
    [I, normI]=BSCBIter(I,Mask,Mask_val, A,GetNorm,t);
    if GetNorm==1
        normsI(1+(switches-1)*(A+B):1+(switches-1)*(A+B)+(A-1))=normI;
    end
    t=t+dt*A;
    
    % TV Inpainting
    %---------------
    disp([num2str(B) ' iterations of TV inpainting']);
    [I,normI]=TVIter(I,Mask,B,1,GetNorm,t); 
    if GetNorm==1
        normsI(switches*A +(switches-1)*B+1:switches*A +(switches-1)*B+B)=normI;
    end
    t=t+dt*B;
end

disp('Applying 10 iterations of TV to entire image');
lambda_val=100;
I=TVIter(I, zeros(size(Mask)), 30,0,0,t,lambda_val);

colormap(gray);
saveas(gcf,strcat('BSCBOutput/', FigNum, 'InpaintedImage.eps'),'epsc');

% Displaying ||I|| vs. time
%----------------------------
if GetNorm==1
    figure;
    plot(0:dt:dt*(length(normsI)-1), normsI,'LineWidth', 1.5);
    title('||Image|| vs. time t', 'FontSize',12);
    xlabel('t');
    ylabel('||Image||')
    saveas(gcf, strcat('BSCBOutput/', FigNum, 'NormI.eps'),'epsc');
end

% Displaying CPU Time
%-----------------------
T=t;                    %final time
CPUtime=toc(tStart)/60;    %total CPU time
disp(' ');
disp(['Results for ' FigNum]);
disp(['Final time t= ' num2str(T)]);
disp(['CPU time= ' num2str(CPUtime) 'mins']);




