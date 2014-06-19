function [I,normI]=TVIter(I0,Mask,MaxIter,MaskOnly,GetNorm, t,lambda_val)
%--------------------------------------------------------------------------
% TVIter performs TV inpainting/denoising using gradient descent
% The number of iterations are determined by the input variable MaxIter
% Inputs:
% I0: initial image, if using for denoising then will be noisy
%                    if using for inpainting, will have masked/unknown regions
% Mask: binary mask for I0 with 1 representing regions to be inpainted
%       If using for denoising, will be a matrix of zeros.
% MaxIter: maximum number of iterations
% MaskOnly: =1 if TV regularization is only applied to masked regions
% GetEnergy: =1 if energy values are computed and stored for each iteration
% lambda_val: optional parameter for specifying lambda for denoising,
%             default value=1e4
% Outputs:
% I: denoised/inpainted image
% normI: ||I|| for each iteration
%--------------------------------------------------------------------------
global dt;

I = I0;
[M,N]=size(I0);
if GetNorm==1; normI=zeros(1,MaxIter); end

% Parameters
%-------------
%  space discretization
h=1.;

% lambda (outside of Mask)
% default value
if nargin==5
    lambda_val=1e4;
end

UpdateI=0;

% Beginning Iterations
%-------------------------
for Iter=1:MaxIter
    for i=2:M-1,
        for j=2:N-1,
            if Mask(i,j)==1
                lambda=0;
                UpdateI=1;
            end
            if Mask(i,j)==0 && MaskOnly==1
                UpdateI=0;
            end
            if Mask(i,j)==0 && MaskOnly==0
                UpdateI=1;
                lambda=lambda_val;
            end
            
            if UpdateI==1
                %-----------computation of coefficients co1,co2,co3,co4---------
                Ix=(I(i+1,j)-I(i,j))/h;
                Iy=(I(i,j+1)-I(i,j-1))/(2*h);
                GradI=sqrt(eps*eps+Ix*Ix+Iy*Iy);
                co1=1./GradI;
                
                Ix=(I(i,j)-I(i-1,j))/h;
                Iy=(I(i-1,j+1)-I(i-1,j-1))/(2*h);
                GradI=sqrt(eps*eps+Ix*Ix+Iy*Iy);
                co2=1./GradI;
                
                Ix=(I(i+1,j)-I(i-1,j))/(2*h);
                Iy=(I(i,j+1)-I(i,j))/h;
                GradI=sqrt(eps*eps+Ix*Ix+Iy*Iy);
                co3=1./GradI;
                
                Ix=(I(i+1,j-1)-I(i-1,j-1))/(2*h);
                Iy=(I(i,j)-I(i,j-1))/h;
                GradI=sqrt(eps*eps+Ix*Ix+Iy*Iy);
                co4=1./GradI;
                
                co=1.+2*lambda*dt +(dt/(h*h))*(co1+co2+co3+co4);
                
                div=co1*I(i+1,j)+co2*I(i-1,j)+co3*I(i,j+1)+co4*I(i,j-1);
                
                I(i,j)=(1./co)*(I(i,j)+2*lambda*dt*I0(i,j)+(dt/(h*h))*div);
            end
            
        end
    end
    
    
    %------------ free boundary conditions -------------------
    for i=2:M-1,
        I(i,1)=I(i,2);
        I(i,N)=I(i,N-1);
    end
    
    for j=2:N-1,
        I(1,j)=I(2,j);
        I(M,j)=I(M-1,j);
    end
    
    I(1,1)=I(2,2);
    I(1,N)=I(2,N-1);
    I(M,1)=I(M-1,2);
    I(M,N)=I(M-1,N-1);
    
    t=t+dt;
    
    imagesc(I);
    %colormap(gray);
    hold on 
    set(gca,'xtick',[],'ytick',[]);
    xlabel(['t= ' num2str(t)],'FontSize',12);
    title('Inpainted Image', 'FontSize',12);
    drawnow;
    hold off;
    
    %---------------- end iterations ------------------------------
    if GetNorm==1; normI(Iter)=norm(I,2); end
end
