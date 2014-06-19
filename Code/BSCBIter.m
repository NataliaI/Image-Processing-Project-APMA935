function[I,normI]=BSCBIter(I,Mask,Mask_val,MaxIter,GetNorm,t)
%--------------------------------------------------------------------------
% BSCBIter performs Bertalmio-Sapiro-Caselles-Ballester inpainting.
% The number of iterations are determined by the input variable MaxIter
% Inputs:
% I: initial image (with masked regions usually marked by a different colour (in white)
% Mask: binary mask with 1 representing regions to be inpainted
% MaxIter: maximum number of iterations
% Outputs:
% I: inpainted image I
% normI: ||I|| for each iteration
%--------------------------------------------------------------------------
global dt; 
[M,N]=size(I);
if GetNorm==1; normI=zeros(1,MaxIter); end
for Iter=1:MaxIter
    for i=1:M
        for j=1:N
            
            % only working inside inpainting region Omega
            %---------------------------------------------
            if Mask(i,j)==1
                % skips values of I that are completely inside the mask
                % updating them won't introduce any new information into
                % the inpainting region
               if I(i,j)==Mask_val && I(i+1,j)==Mask_val && I(i-1,j)==Mask_val && I(i,j+1)==Mask_val && I(i,j-1)==Mask_val...
                       && I(i+2,j)==Mask_val && I(i-2,j)==Mask_val && I(i,j-2)==Mask_val && I(i,j+2)==Mask_val
                   break;
               end
                
                % computing isophote direction N/|N|
                %------------------------------------
                
                % Tried two point and three point stencils for forward and
                % backward differences- don't really make a big difference
                
                % Ix
                %-----
                if I(i+1,j)==Mask_val % if the mask value hasn't been changed then we can figure out what side the boundary is located on
                    %2 pt. backward difference
                    Ix=I(i,j)-I(i-1,j);
                    %3 pt. backward difference
                    %Ix=(3*I(i,j)-4*I(i-1,j)+I(i-2,j))/2;
                end
                if I(i-1,j)==Mask_val
                    %2 pt. forward difference
                    Ix=I(i+1,j)-I(i,j);
                    %3 pt. forward difference
                    %Ix=(-I(i+2,j)+4*I(i+1,j)-3*I(i,j))/2;
                end
                if I(i+1,j)~=Mask_val && I(i-1,j)~=Mask_val
                    % centred difference
                    Ix=(I(i+1,j)-I(i-1,j))/2;
                    %5 pt. centred stencil
                    %Ix=(-I(i+2,j)+8*I(i+1,j)-8*I(i-1,j)+I(i-2,j))/12;
                end
                
                % Iy
                %-----
                if I(i,j+1)==Mask_val
                    %2 pt. backward difference
                    Iy=I(i,j)-I(i,j-1);
                    %3 pt. backward difference
                    %Iy=(3*I(i,j)-4*I(i,j-1)+I(i,j-2))/2;
                end
                if I(i,j-1)==Mask_val
                    %2 pt. forward difference
                    Iy=I(i,j+1)-I(i,j);
                    %3 pt. forward difference
                    %Iy=(-I(i,j+2)+4*I(i,j+1)-3*I(i,j))/2;
                end
                if I(i,j+1)~=Mask_val && I(i,j-1)~=Mask_val
                    %centred difference
                    Iy=(I(i,j+1)-I(i,j-1))/2;
                    % 5 pt. centred difference
                    %Iy=(-I(i,j+2)+8*I(i,j+1)-8*I(i,j-1)+I(i,j-2))/12;
                end
                
                normN=sqrt(Ix^2+Iy^2+eps);
                N_norm=[-Iy Ix]/normN;
                
                % computing change in smoothness of image (dL)
                %-----------------------------------------------
                % dL=[L(i+1,j)-L(i-1,j) L(i,j+1)-L(i,j-1)];
                dL=[laplacian(I,i+1,j)-laplacian(I,i-1,j) laplacian(I,i,j+1)-laplacian(I,i,j-1)];
                
                % computing projection of dL onto N/|N|
                %---------------------------------------
                B=dL*N_norm';
                
                % computing |grad I| using slope-limiters
                %------------------------------------------
                Ixf=I(i+1,j)-I(i,j);  % forward difference
                Ixb=I(i,j)-I(i-1,j);  % backward difference
                Iyf=I(i,j+1)-I(i,j);
                Iyb=I(i,j)-I(i,j-1);
                
                if B>0
                    normGradI=sqrt((min(Ixb,0))^2+(max(Ixf,0))^2+(min(Iyb,0))^2+(max(Iyf,0))^2);
                else
                    normGradI=sqrt((max(Ixb,0))^2+(min(Ixf,0))^2+(max(Iyb,0))^2+(min(Iyf,0))^2);
                end
                
                It=B*normGradI;
                
                % updating I
                %--------------
                I(i,j)=I(i,j)+dt*It;
            end
        end
        
    end
    t=t+dt;
    if GetNorm==1; normI(Iter)=norm(I,2); end % norm of I
    
    % plotting
    %-----------
    imagesc(I); hold on 
    set(gca, 'xtick', [], 'ytick', []);
    xlabel(['t=' num2str(t)], 'FontSize',12);
    title('Inpainted Image', 'FontSize', 12);
    drawnow;
    hold off;
end





