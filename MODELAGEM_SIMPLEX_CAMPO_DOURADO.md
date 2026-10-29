# Modelagem por Programação Linear (Simplex) da dieta de custo mínimo — Faz. Campo Dourado (GO)

**Autor:** Igor Coutinho Ferreira Nobre
**Disciplina:** Pesquisa Operacional — Pós-graduação IFG
**Data:** 25/10/2026

Problema: formular, ao menor custo, a dieta diária de uma vaca do lote de lactação da Fazenda Campo
Dourado, escolhendo quanto de cada alimento disponível entra na ração, de modo que a dieta atenda às
exigências nutricionais do lote (NASEM 2021) e aos limites de formulação. Os dados (alimentos, preços,
composição e exigências) são os da fazenda `GO_campo_dourado/`.

---

## 1. Quais são as variáveis de decisão?

x_j = quantidade do alimento j na dieta, em **kg de matéria seca (MS) por vaca por dia**, j = 1..6.

| j | Alimento | Categoria | c_j (US$/kg MS) |
|---|---|---|---|
| 1 | Silagem de milho | Forragem | 0,16 |
| 2 | Silagem de leguminosa | Forragem | 0,18 |
| 3 | Milho grão moído | Concentrado | 0,19 |
| 4 | Farelo de soja 48 % | Concentrado | 0,41 |
| 5 | Caroço de algodão | Concentrado | 0,36 |
| 6 | Núcleo mineral | Mineral | 0,80 |

---

## 2. Qual é o objetivo?

Minimizar o custo diário da dieta de uma vaca em lactação, atendendo às exigências de energia,
proteína, cálcio, fósforo, fibra, gordura e consumo de matéria seca do lote.

---

## 3. Qual é a função objetivo?

```
min Z = Σ_j c_j x_j = 0,16 x1 + 0,18 x2 + 0,19 x3 + 0,41 x4 + 0,36 x5 + 0,80 x6     [US$/vaca/dia]
```

---

## 4. Quais são as restrições?

### 4.1 Parâmetros

a_ij = quantidade do nutriente i em 1 kg de MS do alimento j (tabela NASEM de composição dos alimentos):

| Nutriente i (por kg MS) | x1 Sil. milho | x2 Sil. legum. | x3 Milho | x4 F. soja | x5 Car. algodão | x6 Mineral |
|---|---|---|---|---|---|---|
| Energia líquida EL_j, Mcal | 1,75 | 1,54 | 2,06 | 2,36 | 1,87 | 0 |
| Proteína metabolizável PM_j, g | 77,8 | 92,1 | 96,8 | 236,5 | 102,6 | 4,7 |
| Cálcio absorvível Ca_j, g | 0,71 | 3,76 | 0,21 | 2,40 | 1,04 | 146,6 |
| Fósforo absorvível P_j, g | 1,48 | 2,25 | 2,16 | 5,17 | 4,35 | 18,56 |
| FDN_j, % | 40,93 | 43,20 | 9,76 | 11,06 | 50,58 | 0,14 |
| FDN de forragem FDNf_j, % | 40,93 | 43,20 | 0 | 0 | 0 | 0 |
| Extrato etéreo (gordura) EE_j, % | 2,92 | 2,87 | 3,85 | 1,81 | 18,62 | 0 |

b_i = exigência diária do lote de lactação (vaca média do lote: 635 kg de peso vivo, 40 kg de leite/dia; média das 41 formulações
da simulação de 5 anos):

| Exigência | b_i |
|---|---|
| Energia líquida | 45,5 Mcal |
| Proteína metabolizável | 2 575 g (máximo permitido: 1,5 × = 3 863 g) |
| Cálcio absorvível | 65,2 g |
| Fósforo absorvível | 61,7 g |
| Consumo de MS estimado (CMS) | 24,4 kg (faixa permitida: 70 % a 130 %) |
| FDN total | entre 25 % e 45 % da MS |
| FDN de forragem | ≥ 15 % da MS |
| Gordura | ≤ 7 % da MS |

### 4.2 Restrições nutricionais (nutriente entregue versus exigência)

```
R1  energia líquida        1,75 x1 + 1,54 x2 + 2,06 x3 + 2,36 x4 + 1,87 x5 + 0 x6        ≥ 45,5
R2  proteína metab. mín.   77,8 x1 + 92,1 x2 + 96,8 x3 + 236,5 x4 + 102,6 x5 + 4,7 x6    ≥ 2 575
R3  proteína metab. máx.   77,8 x1 + 92,1 x2 + 96,8 x3 + 236,5 x4 + 102,6 x5 + 4,7 x6    ≤ 3 863
R4  cálcio                 0,71 x1 + 3,76 x2 + 0,21 x3 + 2,40 x4 + 1,04 x5 + 146,6 x6    ≥ 65,2
R5  fósforo                1,48 x1 + 2,25 x2 + 2,16 x3 + 5,17 x4 + 4,35 x5 + 18,56 x6    ≥ 61,7
```

### 4.3 Restrições em porcentagem da dieta

Fibra e gordura são exigidas como fração da MS total, Σ a_ij x_j / Σ x_j. Como Σ x_j > 0, a razão
é multiplicada pelo consumo e a restrição fica linear:

```
Σ_j a_ij x_j / Σ_j x_j ≥ k   ⇔   Σ_j (a_ij − k) x_j ≥ 0
Σ_j a_ij x_j / Σ_j x_j ≤ k   ⇔   Σ_j (k − a_ij) x_j ≥ 0
```

```
R6  FDN ≥ 25 % da MS          15,93 x1 + 18,20 x2 − 15,24 x3 − 13,94 x4 + 25,58 x5 − 24,86 x6  ≥ 0
R7  FDN ≤ 45 % da MS           4,07 x1 +  1,80 x2 + 35,24 x3 + 33,94 x4 −  5,58 x5 + 44,86 x6  ≥ 0
R8  FDN de forragem ≥ 15 %    25,93 x1 + 28,20 x2 − 15 x3    − 15 x4    − 15 x5    − 15 x6     ≥ 0
R9  gordura ≤ 7 % da MS        4,08 x1 +  4,13 x2 +  3,15 x3 +  5,19 x4 − 11,62 x5 +  7,00 x6  ≥ 0
```

### 4.4 Consumo de matéria seca

```
R10 CMS mínimo    x1 + x2 + x3 + x4 + x5 + x6  ≥  0,70 × 24,4 = 17,0
R11 CMS máximo    x1 + x2 + x3 + x4 + x5 + x6  ≤  1,30 × 24,4 = 31,7
```

### 4.5 Limites por alimento e não negatividade

```
L_j ≤ x_j ≤ U_j,   x_j ≥ 0,   j = 1..6
```

- **Otimização livre:** L_j = 0 e U_j = 100 kg (sem limite prático).
- **Receita do usuário (modo em uso na fazenda):** cada alimento fica a ±10 % da participação ρ_j
  definida pelo produtor (30 / 18 / 24 / 16 / 9 / 3 % da MS), aplicada sobre a meta de consumo
  1,1 × CMS = 26,84 kg:

  ```
  L_j = 0,9 × ρ_j × 26,84        U_j = 1,1 × ρ_j × 26,84
  ```

  | | x1 | x2 | x3 | x4 | x5 | x6 |
  |---|---|---|---|---|---|---|
  | ρ_j | 0,30 | 0,18 | 0,24 | 0,16 | 0,09 | 0,03 |
  | L_j (kg MS) | 7,25 | 4,35 | 5,80 | 3,87 | 2,17 | 0,72 |
  | U_j (kg MS) | 8,86 | 5,31 | 7,09 | 4,72 | 2,66 | 0,89 |

---

## 5. Modelo completo

```
min  Z = 0,16 x1 + 0,18 x2 + 0,19 x3 + 0,41 x4 + 0,36 x5 + 0,80 x6

s.a. Σ_j EL_j x_j              ≥ 45,5        (R1  energia)
     Σ_j PM_j x_j              ≥ 2 575       (R2  proteína mín.)
     Σ_j PM_j x_j              ≤ 3 863       (R3  proteína máx.)
     Σ_j Ca_j x_j              ≥ 65,2        (R4  cálcio)
     Σ_j P_j  x_j              ≥ 61,7        (R5  fósforo)
     Σ_j (FDN_j − 25) x_j      ≥ 0           (R6  FDN mín.)
     Σ_j (45 − FDN_j) x_j      ≥ 0           (R7  FDN máx.)
     Σ_j (FDNf_j − 15) x_j     ≥ 0           (R8  FDN forragem)
     Σ_j (7 − EE_j) x_j        ≥ 0           (R9  gordura)
     Σ_j x_j                   ≥ 17,0        (R10 CMS mín.)
     Σ_j x_j                   ≤ 31,7        (R11 CMS máx.)
     L_j ≤ x_j ≤ U_j,  x_j ≥ 0               (limites por alimento)
```

6 variáveis de decisão e 11 restrições funcionais, mais os limites por alimento.

| Restrição | x1 | x2 | x3 | x4 | x5 | x6 | sinal | b |
|---|---|---|---|---|---|---|---|---|
| R1 energia | 1,75 | 1,54 | 2,06 | 2,36 | 1,87 | 0 | ≥ | 45,5 |
| R2 PM mín. | 77,8 | 92,1 | 96,8 | 236,5 | 102,6 | 4,7 | ≥ | 2 575 |
| R3 PM máx. | 77,8 | 92,1 | 96,8 | 236,5 | 102,6 | 4,7 | ≤ | 3 863 |
| R4 cálcio | 0,71 | 3,76 | 0,21 | 2,40 | 1,04 | 146,6 | ≥ | 65,2 |
| R5 fósforo | 1,48 | 2,25 | 2,16 | 5,17 | 4,35 | 18,56 | ≥ | 61,7 |
| R6 FDN mín. | 15,93 | 18,20 | −15,24 | −13,94 | 25,58 | −24,86 | ≥ | 0 |
| R7 FDN máx. | 4,07 | 1,80 | 35,24 | 33,94 | −5,58 | 44,86 | ≥ | 0 |
| R8 FDN forr. | 25,93 | 28,20 | −15 | −15 | −15 | −15 | ≥ | 0 |
| R9 gordura | 4,08 | 4,13 | 3,15 | 5,19 | −11,62 | 7,00 | ≥ | 0 |
| R10 CMS mín. | 1 | 1 | 1 | 1 | 1 | 1 | ≥ | 17,0 |
| R11 CMS máx. | 1 | 1 | 1 | 1 | 1 | 1 | ≤ | 31,7 |
| **custo c_j** | 0,16 | 0,18 | 0,19 | 0,41 | 0,36 | 0,80 | min | |

Observação sobre os coeficientes de energia e proteína: no NASEM eles dependem da dieta inteira
(consumo total, peso vivo, proporção de concentrado). Os valores da tabela são os coeficientes
lineares válidos em torno da dieta praticada na fazenda (linearização em torno do ponto de operação);
quando a solução se afasta desse ponto, os coeficientes são recalculados e o PL é resolvido de novo.
