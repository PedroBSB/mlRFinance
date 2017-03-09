mlRFinance
---
Um pacote em R para geração de portifólios baseado em aprendizado de máquinas. 


Instalação
===
######Automatica:
1. No R instale o devtools ```install.packages("devtools")```
2. Após finalizada a instalação que pode demorar alguns minutos instale o seguinte pacote do github ```devtools::install_github("hadley/devtools")```
3. Finalizado o processo execute ```devtools::install_github("PedroBSB/mlRFinance", auth_token="0c468b0e44dd1753f80ba8fe7e73d5187b98563a")```

######Manual:
1. No Windows, baixe e instale <a href="https://cran.r-project.org/bin/windows/Rtools/" target="_blank">Rtools</a>. 
No Mac, instale o <a href="https://developer.apple.com/xcode/features/" target="_blank">Xcode command line tools</a>.
No Linux, instale o R development package, normalmente chamado de r-devel ou r-base-dev dependendo da sua distribuição.

2. No Rstudio execute o comando ```install.packages("devtools")``` em seguida ```devtools::install.packages("hadley/devtools")```

3. Apos isso instale ou certifique-se que os seguintes pacotes estão instalados: ```Rcpp RcppArmadillo RcppParallel RcppEigen RcppProgress foreach doParallel foreach parallel iterators Matrix boot doSNOW```. Basta usar ```install.packages(c("Rcpp","RcppEigen","RcppArmadillo","RcppParallel","RcppProgress","foreach","doParallel","parallel","iterators","Matrix","boot","doSNOW"),dependencies=TRUE)```

4. Abra o pacote através da opção disponível no Rstudio, após o carregamento dos arquivos do pacote na janela a direta clique em Build e depois no botão 'Build & Reload' ou usando a tecla de atalho Ctrl+Shift+B.

To Do List:
===
- [ ] Converter a função em src\Utils.cpp\nearPDefinite totalmente em RcppEigen.
- [ ] Escrever um Wiki page sobre as medidas de mensuração de erro em src\ErrorMeasure.cpp.
- [ ] Escrever um Wiki page com exemplos das funções (um bom começo é usar os scripts que estão na pasta test).
- [ ] Escrever o manual das funções.
- [ ] Adicionar mais Similarity Matrices na função QPFS.


-----
