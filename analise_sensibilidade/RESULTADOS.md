# Análise de sensibilidade - números para o relatório

Fonte: `output/simplex/GO_campo_dourado/modelo_LAC_COW_dia00045_tent1.csv` e `tableau_LAC_COW_dia00045_tent1.csv`
(formulação do dia 45). Calculado por `sensibilidade.py` (`faixa_custos`, `faixa_rhs`).

## Solução ótima (dia 45)

Z* = 6,1516 US$/vaca/dia, 2 pivôs, 2 iterações da linearização.
Base: x1 (silagem de milho), excessos de R2 a R11, folgas dos 6 tetos.

| | x* (kg MS) | piso L_j | teto U_j | situação |
|---|---|---|---|---|
| x1 silagem de milho | 8,089 | 7,286 | 8,905 | básica |
| x2 silagem de leguminosa | 4,372 | 4,372 | 5,343 | no piso |
| x3 milho grão | 5,829 | 5,829 | 7,124 | no piso |
| x4 farelo de soja | 3,886 | 3,886 | 4,749 | no piso |
| x5 caroço de algodão | 2,186 | 2,186 | 2,672 | no piso |
| x6 núcleo mineral | 0,729 | 0,729 | 0,891 | no piso |

## Preços-sombra e intervalo do lado direito

| restrição | b_i | valor em x* | folga | y_i | b_i mín | b_i máx |
|---|---|---|---|---|---|---|
| R1 energia (Mcal) | 50,981 | 50,981 | 0 | 0,0814 | 49,403 | 52,584 |
| R2 proteína mín (g) | 2727,0 | 2844,6 | 117,6 | 0 | -inf | 2844,6 |
| R3 proteína máx (g) | 4015,8 | 2844,6 | 1171,2 | 0 | 2844,6 | +inf |
| R4 cálcio (g) | 66,6 | 141,8 | 75,2 | 0 | -inf | 141,8 |
| R5 fósforo (g) | 62,8 | 77,6 | 14,8 | 0 | -inf | 77,6 |
| R6 FDN >= 25 % | 0 | 103,2 | 103,2 | 0 | -inf | 103,2 |
| R7 FDN <= 45 % | 0 | 398,6 | 398,6 | 0 | -inf | 398,6 |
| R8 FDN forragem >= 15 % | 0 | 143,6 | 143,6 | 0 | -inf | 143,6 |
| R9 gordura <= 7 % | 0 | 69,3 | 69,3 | 0 | -inf | 69,3 |
| R10 CMS mín (kg) | 17,17 | 25,09 | 7,92 | 0 | -inf | 25,09 |
| R11 CMS máx (kg) | 31,89 | 25,09 | 6,80 | 0 | 25,09 | +inf |

Tetos max_j: todos com folga, y = 0. Pisos: y do piso = custo reduzido do alimento.
Conferência do y da energia: 0,16 / 1,9647 Mcal/kg (silagem) = 0,0814.

Verificação re-otimizando R1:

| b1 | Z previsto = Z* + 0,0814·Δ | Z re-otimizado | base |
|---|---|---|---|
| 49,5 | 6,0310 | 6,0310 | igual |
| 52,0 | 6,2346 | 6,2346 | igual |
| 53,0 | 6,3161 | 6,3187 | muda (x3 entra) |
| 54,0 | 6,3975 | 6,4064 | muda |

## Intervalo de c_j e custos reduzidos

| | c_j | custo reduzido | c_j mín | c_j máx |
|---|---|---|---|---|
| x1 silagem de milho | 0,16 | 0 | 0,00 | 0,1724 |
| x2 silagem de leguminosa | 0,18 | 0,0346 | 0,1454 | +inf |
| x3 milho grão | 0,19 | 0,0137 | 0,1763 | +inf |
| x4 farelo de soja | 0,41 | 0,1983 | 0,2117 | +inf |
| x5 caroço de algodão | 0,36 | 0,1831 | 0,1769 | +inf |
| x6 núcleo mineral | 0,80 | 0,8206 | 0,00 | +inf |

Nota: mineral tem EL = -0,2532 Mcal/kg no modelo linearizado, por isso custo reduzido > preço.
EL por kg MS (linha R1): x1 1,9647, x2 1,7848, x3 2,1651, x4 2,5995, x5 2,1725, x6 -0,2532.

Verificação re-otimizando c_j:

| mudança | Z previsto | Z re-otimizado | dieta |
|---|---|---|---|
| x1 a 0,10 | 5,6663 | 5,6663 | igual |
| x1 a 0,17 | 6,2325 | 6,2325 | igual |
| x1 a 0,18 | 6,3134 | 6,3073 | muda: x1 7,286, x3 6,558 |
| x3 a 0,175 | - | 6,0632 | muda: x1 7,286, x3 6,558 |
| x3 a 0,30 | 6,7928 | 6,7928 | igual |

## Cenário milho +60 % (0,19 -> 0,30), dia 45

| | Z (US$/vaca/dia) | dieta (kg MS) |
|---|---|---|
| antes | 6,1516 | atual |
| manter dieta, preço novo | 6,7928 | atual |
| re-otimizar, receita ±10 % | 6,7928 | atual (x3 já no piso) |
| re-otimizar, receita ±30 % | 6,2723 | x3 4,534, x1 10,524, x5 1,700, x2 5,011, x4 3,076, x6 0,567 (7 pivôs) |

## Dia 1800 (caso degenerado)

Z* = 5,9519, 0 pivôs. Seis alimentos no piso, nenhuma restrição funcional ativa (folga da energia 0,002 Mcal),
todos os y = 0, custo reduzido de cada x_j = c_j, faixa de c_j = [0, +inf) para todos.
x* = 5,76 / 7,20 / 2,16 / 4,32 / 3,84 / 0,72 (x3, x1, x5, x2, x4, x6).

## Preço-sombra da energia nas 41 formulações

| dias | restrição ativa | y |
|---|---|---|
| 0 | proteína mín | 0,0017 US$/g |
| 45 a 585 | energia | 0,0811 a 0,0818 US$/Mcal |
| 630 | energia + proteína mín (5 pivôs) | 0,0260 US$/Mcal e 0,0014 US$/g |
| 675 | energia | 0,0884 US$/Mcal (máximo) |
| 720 a 1755 | energia | 0,0809 a 0,0819 US$/Mcal |
| 1800 | nenhuma | 0 |
