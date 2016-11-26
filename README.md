mlRFinance
=====
Um pacote em R para geração de portifólios baseado em aprendizado de máquinas. 


Instalação
-----

1.No Windows, baixe e instale Rtools.
No Mac, instale o Xcode command line tools.
No Linux, instale o R development package, normalmente chamado de r-devel ou r-base-dev dependendo da sua distribuição.

2.No Rstudio execute o comando install.packages("hadley/devtools")

3.Abra o pacote através da opção disponivel no Rstudio, após o carregamento dos arquivos do pacote na janela a direta clique em Build.


É preciso fazer a importação do LAPACK and BLAS
1) Rode a biblioteca library(RcppArmadillo)
2) Execute RcppArmadillo.package.skeleton("Pacote")
3) Copie os arquivos Makevars para a pasta src do seu pacote
#Olhar aqui:: http://stackoverflow.com/questions/28754573/c-compiling-error-while-compiling-r-package-on-winbuild

#Para o RcppParallel:
https://rcppcore.github.io/RcppParallel/#r_packages
