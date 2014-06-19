% BSCB Inpainting Examples
%---------------------------
% Performs BSCB inpainting using the algorithm outlined in Bertalmio et
% al.'s paper.
% Inpainting regions in I are given a value of 1.3 (outside the colour range)
% so that the algorithm can keep track of the values that have been updated
% in the inpainting region.
% Each cell can be run one at a time to get individual results for each
% figure, or the file can be run all at once as well.
% Figure numbers correspond to project report.

% Calls the functions BSCBInpaint, BSCBInpaintColour, BSCBIter, TVIter,
% laplacian.

clear; close all; clc;
% For all examples
%-------------------
global dt;
dt=0.1;

%% Figure 1: 
%----------------------------------------

I=imread('Images/chevron.jpg');
I=im2double(I);
load Images/chevronMaskMiddle
I(Mask==1)=1.3;

[TFig6, CPUFig6]=BSCBInpaint(I,Mask,1000,1,'Fig6');

%% Figure 2: Mama
%--------------------
I=imread('Images/Old.jpg');
I=im2double(rgb2gray(I));
load Images/OldMask
I(Mask==1)=1.3;

[TFig7, CPUtimeFig7]=BSCBInpaint(I,Mask,550,1,'Fig7');

%% Figure 3: 
I=imread('Images/text.jpg');
I=im2double(rgb2gray(I));
load Images/textMask
I(Mask==1)=1.3;

[TFig8,CPUtimeFig8]=BSCBInpaint(I,Mask,500,1, 'Fig8');

%% Figure 4: COLOUR
I=imread('Images/parachute.jpg');
I=im2double(I);
load Images/parachuteMask

[TFig9,CPUtimeFig9]=BSCBInpaintColour(I,Mask,300, 'Fig9');

%% Figure 5: Einstein
I=imread('Images/einstein.jpg');
I=im2double(I);
load Images/einsteinMask

[TFig10,CPUtimeFig10]=BSCBInpaintColour(I,Mask,100,'Fig10');
