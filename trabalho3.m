1;
function V = SolvePotentialWithIteration(S, kmax) # A decomposição QR Gram-Schmidt modificada
  xSize = size(S,1);
  ySize = size(S,2);
  %VMenos1 = zeros(xSize,ySize);
  for i = 1:xSize
    for j =1:ySize
      V(i,j) = S(i,j).valor;
    endfor
  endfor
  %V = VMenos1;
  k = 0;
  while k < kmax
    for i = 1:xSize
      for j =1:ySize
        if S(i, j).inicial == 0
          soma = 0;

          if i - 1 > 0
            soma += V(i - 1, j);
          endif
          if i + 1 < xSize
            soma += V(i + 1, j);
          endif
          if j - 1 > 0
            soma += V(i, j - 1);
          endif
          if j + 1 < ySize
            soma += V(i, j + 1);
          endif

          V(i,j) = soma/4;
        endif
      endfor
    endfor
  %VMenos1 = V;
  k++
  endwhile

endfunction
function S = ImagemParaProblema(I, significadoDasCores) # A decomposição QR Gram-Schmidt modificada
  xSize = size(I,1);
  ySize = size(I,2);
  %S = zeros(xSize,ySize);
  for i = 1:xSize
    for j =1:ySize
      corID = ['c'  num2str(I(i,j,1))  '-'  num2str(I(i,j,2))  '-'  num2str(I(i,j,3))];
      if isfield(significadoDasCores,corID)
        valorInicial = getfield(significadoDasCores,corID);
        S(i,j) = struct("inicial",true, "valor",valorInicial);
      else
        S(i,j) = struct("inicial",false, "valor",0);
      endif
    endfor
  endfor
endfunction
# Descomentar essa parte caso queira rodar o programa de forma limpa
clear
clc
significadoDasCores = struct("c1-0-1",20, "c0-0-1",0);%, "c0-1-0",40, "c1-0-0",10);
Imagem = imread('problema3.png');
A = ImagemParaProblema(Imagem,significadoDasCores);
%A = [struct("inicial",true, "valor",40) struct("inicial",true, "valor",40) struct("inicial",true, "valor",40) struct("inicial",true, "valor",40);
%     struct("inicial",true, "valor",20) struct("inicial",false, "valor",0) struct("inicial",false, "valor",0) struct("inicial",true, "valor",0);
%     struct("inicial",true, "valor",20) struct("inicial",false, "valor",0) struct("inicial",false, "valor",0) struct("inicial",true, "valor",0);
%     struct("inicial",true, "valor",10) struct("inicial",true, "valor",10) struct("inicial",true, "valor",10) struct("inicial",true, "valor",10)
%    ];
V = SolvePotentialWithIteration(A,300);

xSize = size(V,1);
ySize = size(V,2);
tx = 1:1:xSize;
ty = 1:1:ySize;
%[xx, yy] = meshgrid (tx, ty);

mesh (tx, ty, V');
xlabel ("x");
ylabel ("y");
zlabel ("V");
title ("3-D Potential plot");

#-----------------------------------------------
