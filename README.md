mlRFinance
---
Um pacote em R para gera??o de portif?lios baseado em aprendizado de m?quinas. 


Instala??o
===
######Automatica:
1. No R instale o devtools ```install.packages("devtools")```
2. Ap?s finalizada a instala??o que pode demorar alguns minutos instale o seguinte pacote do github ```devtools::install_github("hadley/devtools")```
3. Finalizado o processo execute ```devtools::install_github("PedroBSB/mlRFinance", auth_token=*chave*)```

######Manual:
1. No Windows, baixe e instale Rtools.
No Mac, instale o Xcode command line tools.
No Linux, instale o R development package, normalmente chamado de r-devel ou r-base-dev dependendo da sua distribui??o.

2. No Rstudio execute o comando ```install.packages("devtools")``` em seguida ```devtools::install.packages("hadley/devtools")```

3. Apos isso instale ou certifique-se que os seguintes pacotes est?o instalados: ```Rcpp RcppArmadillo RcppParallel RcppEigen RcppProgress foreach doParallel```

4. Abra o pacote atrav?s da op??o disponivel no Rstudio, ap?s o carregamento dos arquivos do pacote na janela a direta clique em Build e depois no bot?o 'Build & Reload'.


-----

É preciso fazer a importação do LAPACK and BLAS:

1. Rode a biblioteca library(RcppArmadillo)

2. Execute RcppArmadillo.package.skeleton("Pacote")

3. Copie os arquivos Makevars para a pasta src do seu pacote


**Olhar aqui**: http://stackoverflow.com/questions/28754573/c-compiling-error-while-compiling-r-package-on-winbuild
_____
**Para o RcppParallel**: https://rcppcore.github.io/RcppParallel/#r_packages
