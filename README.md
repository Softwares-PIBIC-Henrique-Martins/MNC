# MNC
Trata-se de um software capaz de realizar ajuste de curvas experimentais obtidas por uma balança analítica para a massa aparente de amostras de nanopartículas magnéticas em função do tempo. O ajuste pode ser feito a partir de três métodos: Método dos Mínimos Quadrados, Método de Monte Carlo e Algoritmo Genético. A partir desse ajuste, o software é capaz de obter estimativas para a constante de anisotropia magnética efetiva (*Kef*) da amostra analisada.

## Sobre os métodos de ajuste
Cada um desses três métodos, para ser executado, precisa de uma estimativa inicial de *Kef* (***Estimativa de Kef***).
O método de ajuste por Monte Carlo realiza ***Nº de ensaios*** ensaios em um intervalo definido por [***Estimativa de Kef - Raio, Estimativa de Kef + Raio***]. e retorna o valor para *Kef* correspondente ao ensaio em que se teve o menor valor para o qui-quadrado.

O método de ajuste por Algoritmo Genético também realiza ***Nº de ensaios*** ensaios em um intervalo definido por [***Estimativa de Kef - Raio, Estimativa de Kef + Raio***]., no entanto, a quantidade de ensaios a serem realizados depende da entrada ***Limite***. O método busca fazer ensaios de maneira indefinida, mas se a ele fizer uma quantidade de ensaios maior que o número ***Limite*** sem conseguir achar um valor melhor para *Kef* (isso é, que corresponde a um qui-quadrado menor), então a função que faz o ajuste para e então já teremos nosso valor estimado para o *Kef*.

## inputs
- ***Dados:*** os dados das curvas experimentais podem ser digitados diretamente na caixa de texto contida no frame **Dados** ou então podem ser inseridos a partir de uma pesquisa no diretório ao acessar *Arquivo > Abrir*. O arquivo deve estar em formato *xlsx*, com uma coluna intitulada "t" para os dados do tempo e uma outra coluna intitulada "m", para os dados da massa aparente. Após os dados serem inseridos na caixa de texto, podem ser editados.

- ***Parâmetros:*** os parâmetros variáveis do experimento se encontram no frame **Parâmetros**, onde o usuário deve especificar o valor de cada parâmetro para que o ajuste possa ser feito.

- ***Tipo de ajuste:*** o usuário deve clicar nas checkboxes correspondentes aos métodos de ajustes que ele quer executar. Essas checkboxes se encontram no frame **Plotagem**.

- ***Inputs dos ajustes:*** o usuário deve digitar o valor das entradas das funções de ajuste, como a ***Estimativa de *Kef***, ***Raio***, ***Nº de ensaios*** e ***Limite***. Esses valores devem ser digitados nas entrys presentes no frame **Plotagem**.

## outputs
Os retornos dados pelo software após todas as inputs serem informadas adequadamente e após o botão **Plotar** ser clicado são:
- ***Kef***: o valor obtido pelo ajuste para *Kef*. Esse valor aparece no frame **Saída**
- ***Erro***: é o valor obtido para o erro associado ao valor obtido para *Kef*. Esse retorno só é associado ao método dos Mínimos Quadrados e aparece no frame **Saída**.
- ***Quiquadrado***: é o valor do quiquadrado associado ao valor obtido para *Kef* para cada ajuste. Esse valor aparece no frame **Saída**
- ***Gráfico***: há a plotagem da curva experimental e da curva de ajuste.
