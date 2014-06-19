%  TV Inpainting Examples
%---------------------------
% Performs TV Inpainting using gradient descent, with semi-implicit Gauss-
% Seidel iteration. 
% All masks were created beforehand.
% Regions to be inpainted in the image,I, are set to whatever colour value desired
% Most of the time setting it equal to random values in [0,1] produces the
% fastest results. 
% Each cell can be run one at a time to get individual results for each
% figure, or the file can be run all at once as well.
% Figure 4 takes about 15mins to run. 

% Calls the functions TVInpaint, TVtol.

clear; close all; clc;
% For all examples: 
%-------------------------
global dt;
dt=0.1;

%% Figure1: Bar Image with Thin Mask and Thick Mask
%-------------------------------------------
    I=imread('Images/bar.jpg');
    I=im2double(I);

    % Thin Mask 1a)
    %------------
    load Images/barMaskThin
    I(Mask==1)=rand(nnz(Mask==1),1);
    tol=1e-3;
    TVInpaint(I,Mask,tol,'Fig1a');
    
    % Thick Mask 1b)
    %------------
    load Images/barMaskWide
    I(Mask==1)=rand(nnz(Mask==1),1);
    TVInpaint(I,Mask,tol,'Fig1b');
    
%% Figure 2: Chevron Image with 2 Masks
%----------------------------------
    I=imread('Images/chevron.jpg');
    I=im2double(I);
    tol=1e-4;
    load Images/chevronMaskBottom
    I(Mask==1)=rand(nnz(Mask==1),1);
    TVInpaint(I,Mask,tol,'Fig2a');

    I=imread('Images/chevron.jpg');
    I=im2double(I);
    tol=1e-4;
    load Images/chevronMaskMiddle
    I(Mask==1)=rand(nnz(Mask==1),1);
    TVInpaint(I,Mask,tol,'Fig2b');
    
%% Figure 3: Storm
%--------------------
    I=imread('Images/Storm.jpg');
    I=im2double(rgb2gray(I));
    tol=1e-1;
    load Images/StormMask
    I(Mask==1)=1;
    figure;imagesc(I); colormap(gray);
    TVInpaint(I,Mask,tol,'Fig3');

