function[]=TVInpaint(I,Mask,tol,FigNum)
%--------------------------------------------------------------------------
% TVInpaint calls the function TVtol which performs the actual inpainting
% algorithm
% The initial, final, and energy figures are all plotted and output to the
% file 'Outputs'.
% Inputs: I:      initial image
%         Mask:   binary mask of region to be inpainted
%         tol:    tolerance level of Energy
%         FigNum: figure number for ouput (string format)
%--------------------------------------------------------------------------
global dt

% Initial Image
%----------------
figure;
imagesc(I); axis off; colormap(gray);
title('Initial Image', 'FontSize', 12);
saveas(gcf, strcat('TVOutput/',FigNum,'InitialImage.eps'),'epsc');

% Inpainting
%-------------
tStart=tic;
figure;
[I,t,Energy]=TVtol(I,Mask,tol,1); %TVtol(I0, Mask, tol_val, MaskOnly)
CPUTime=toc(tStart);
saveas(gca, strcat('TVOutput/',FigNum,'InpaintedImage.eps'),'epsc');

% Energy vs. Iterations
%-------------------------
figure;
plot(0:dt:dt*(length(Energy)-1),Energy,'LineWidth',1.5);
title(['Energy vs. Iterations with tol= ' num2str(tol)],'FontSize',12);
ylabel('Energy');
xlabel('t');
saveas(gcf,strcat('TVOutput/',FigNum,'EnergyPlot.eps'),'epsc');

% Displaying CPU Time
%----------------------
disp(' ');
disp(['Results for ' FigNum]);
disp('-----------------------------');
disp(['Final time t= ' num2str(t)]); 
disp(['CPU Time: ' num2str(CPUTime/60) ' mins']);