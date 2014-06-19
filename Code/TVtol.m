function [I,t,Energy]=TVtol(I0,Mask,tol_val, MaskOnly, lambda_val)
%--------------------------------------------------------------------------
% TVtol performs TV inpainting/denoising using gradient descent
% The number of iterations are determined by setting a tolerance level on
% the energy.
% Inputs:
% I0: initial image, if using for denoising then will be noisy
%                    if using for inpainting, will have masked/unknown regions
% Mask: binary mask for I0 with 1 representing regions to be inpainted
%       If using for denoising, will be a matrix of zeros.
% tol_val: tolerance level for energy
% MaskOnly: =1 if TV regularization is only applied to masked regions
% lambda_val: optional parameter for specifying lambda for denoising,
%             default value=1e4
% Outputs:
% I: denoised/inpainted image
% t: final time t
% Energy: value of energy at each iteration
%--------------------------------------------------------------------------
global dt

I = I0;
[M,N]=size(I0);

% Parameters
%--------------
%  space discretization
h=1.;

% lambda (outside of Mask)
% default
if nargin==4
    lambda_val=1e4;
end

UpdateI=0;

% Beginning Iterations
%-----------------------
t=0;
tol=10000;
Iter=1;
while tol>tol_val
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
    
    %------------- plotting -------------------------------------
    imagesc(I);
    colormap(gray);
    hold on 
    set(gca,'xtick',[],'ytick',[]);
    xlabel(['t= ' num2str(t)],'FontSize',12);
    title('Inpainted Image', 'FontSize',12);
    drawnow;
    hold off;
    
    %---------------- end iterations ------------------------------
    
    % Compute the discrete energy at each iteration
    %-----------------------------------------------
    en=0.0;
    for i=2:M-1,
        for j=2:N-1,
            Ix=(I(i+1,j)-I(i,j))/h;
            Iy=(I(i,j+1)-I(i,j))/h;
            fid=(I0(i,j)-I(i,j))*(I0(i,j)-I(i,j));
            en=en+sqrt(eps*eps+Ix*Ix+Iy*Iy)+lambda*fid;
        end
    end
    
    Energy(Iter)=en;
    
    if Iter>1
        tol=abs(Energy(Iter)-Energy(Iter-1));
    end
    
    t=t+dt;
    Iter=Iter+1;
end
