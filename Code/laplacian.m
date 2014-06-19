function [L]=laplacian(I,i,j)

    L=I(i-1,j)+I(i+1,j)+I(i,j-1)+I(i,j+1)-4*I(i,j);
end