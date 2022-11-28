1;
function V = SolvePotentialWithLinearSystem(S)
  xSize = size(S,1);
  ySize = size(S,2);
  numberOfVariables = xSize * ySize;
  equationCoeficients = zeros(numberOfVariables, numberOfVariables);
  equationValues = zeros(numberOfVariables,1);
  V = zeros(xSize,ySize);

  for j = 1:ySize
    for i =1:xSize
      #Equação da EDP:
      #4 * V(posInEquation) - V(posInEquationOfUpValue) - V(posInEquationOfDownValue) - V(posInEquationOfLeftValue) - V(posInEquationOfRightValue) = 0
      #Se posInEquation é inicial, não mudar a matriz
      #Se posInEquationOfUpValue, posInEquationOfDownValue, posInEquationOfLeftValue ou posInEquationOfRightValue forem iniciais
      #no lugar de fazer com que equationValues(posInEquation) = 0, será igual a soma dos inciais
      #e se eles sairem dos limites da matriz, falar que eles são iniciais e iguais a zero, ou seja, ignorar eles
      if S(i,j).inicial == 0
        posInEquation =  i + (j - 1) * ySize;
        posInEquationOfUpValue =  i + (j - 2) * ySize;
        posInEquationOfDownValue =  i + j * ySize;
        posInEquationOfLeftValue =  (i - 1) + (j - 1) * ySize;
        posInEquationOfRightValue =  (i + 1) + (j - 1) * ySize;

        lineValue = 0;
        equationCoeficients(posInEquation, posInEquation) = 4;
        if (i  <= xSize && i  >= 1) && (j - 1 <= ySize && j - 1 >= 1)
          if S(i,j - 1).inicial == 0
            equationCoeficients(posInEquation, posInEquationOfUpValue) = -1;
          else
            lineValue += S(i,j - 1).valor;
          endif
        endif

        if (i <= xSize && i >= 1) && (j + 1 <= ySize && j + 1 >= 1)
          if S(i,j + 1).inicial == 0
            equationCoeficients(posInEquation, posInEquationOfDownValue) = -1;
          else
            lineValue += S(i,j + 1).valor;
          endif
        endif

        if (i - 1 <= xSize && i - 1 >= 1) && (j <= ySize && j >= 1)
          if S(i - 1,j).inicial == 0
            equationCoeficients(posInEquation, posInEquationOfLeftValue) = -1;
          else
            lineValue += S(i - 1,j).valor;
          endif
        endif

        if (i + 1 <= xSize && i + 1 >= 1) && (j <= ySize && j >= 1)
          if S(i + 1,j).inicial == 0
            equationCoeficients(posInEquation, posInEquationOfRightValue) = -1;
          else
            lineValue += S(i + 1,j).valor;
          endif
        endif
      equationValues(posInEquation) = lineValue;
      endif
    endfor
  endfor
  #equationValues
  #equationCoeficients
  result = pinv(equationCoeficients) * equationValues;
  #result
  for j = 1:ySize
    for i =1:xSize
      if S(i,j).inicial == 0
        posInEquation =  i + (j - 1) * ySize;
        #result(posInEquation)
        V(i,j) = result(posInEquation);
      else
        V(i,j) = S(i,j).valor;
      endif
    endfor
  endfor
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
Imagem = imread('problema.png');
A = ImagemParaProblema(Imagem,significadoDasCores);
#A = [struct("inicial",true, "valor",40) struct("inicial",true, "valor",40) struct("inicial",true, "valor",40) struct("inicial",true, "valor",40);
#     struct("inicial",true, "valor",20) struct("inicial",false, "valor",0) struct("inicial",false, "valor",0) struct("inicial",true, "valor",0);
#     struct("inicial",true, "valor",20) struct("inicial",false, "valor",0) struct("inicial",false, "valor",0) struct("inicial",true, "valor",0);
#     struct("inicial",true, "valor",10) struct("inicial",true, "valor",10) struct("inicial",true, "valor",10) struct("inicial",true, "valor",10)
#    ];
V = SolvePotentialWithLinearSystem(A);

xSize = size(V,1);
ySize = size(V,2);
tx = 1:1:xSize;
ty = 1:1:ySize;
[xx, yy] = meshgrid (tx, ty);
mesh (tx, ty, V');

#contour(tx, ty, V', 0:20);

xlabel ("x");
ylabel ("y");
zlabel ("V");
title ("3-D Potential plot");

#-----------------------------------------------
