mlRFinance
---
Um pacote em R para geração de portifólios baseado em aprendizado de máquinas. 


Instalação
===
######Automatica:
1. No R instale o devtools ```install.packages("devtools")```
2. Após finalizada a instalação que pode demorar alguns minutos instale o seguinte pacote do github ```devtools::install_github("hadley/devtools")```
3. Finalizado o processo execute ```devtools::install_github("PedroBSB/mlRFinance", auth_token=*chave*)```

######Manual:
1. No Windows, baixe e instale Rtools.
No Mac, instale o Xcode command line tools.
No Linux, instale o R development package, normalmente chamado de r-devel ou r-base-dev dependendo da sua distribuição.

2. No Rstudio execute o comando ```install.packages("devtools")``` em seguida ```devtools::install.packages("hadley/devtools")```

3. Apos isso instale ou certifique-se que os seguintes pacotes estão instalados: ```Rcpp RcppArmadillo RcppParallel RcppEigen RcppProgress foreach doParallel foreach parallel iterators Matrix boot```. Basta usar ```install.packages("Rcpp","RcppEigen","RcppArmadillo","RcppParallel","RcppProgress","foreach","doParallel","parallel","iterators","Matrix","boot",dependencies=TRUE)```

4. Abra o pacote através da opção disponível no Rstudio, após o carregamento dos arquivos do pacote na janela a direta clique em Build e depois no botão 'Build & Reload' ou usando a tecla de atalho Ctrl+Shift+B.


-----
