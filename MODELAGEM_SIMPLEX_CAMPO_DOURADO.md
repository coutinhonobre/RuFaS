# Modelagem por Programação Linear (Simplex) da dieta de custo mínimo — Faz. Campo Dourado (GO)

**Autor:** Igor Coutinho Ferreira Nobre
**Disciplina:** Pesquisa Operacional — Pós-graduação IFG
**Data:** 25/10/2026

Projeto final de Pesquisa Operacional: *"Desenvolvimento de módulo de otimização da dieta
integrada ao RuFaS com avaliação de efeito econômico, e também análise de sensibilidade."*

Este documento cobre **só a primeira etapa: a modelagem**, como Programação Linear (PL),
para ser resolvida pelo **Simplex**. Não há código aqui. A estrutura segue o quadro do professor
(`material_complementar/DescricaoDamodelagem.jpeg`), que divide o problema em três blocos:

```
 (1) Estágio do animal  →  (2) Requisitos nutricionais  →  (3) Otimização
     BW, lactação, DEL,        NASEM (2021) / NRC (2001):      ingredientes, kg, $/ton
     produção de leite          NEL ≥ …, FDN ≤ …, …            → $/vaca/dia
```

O RuFaS faz exatamente esse caminho, por lote (pen) de animais, a cada intervalo de formulação.
Fontes de dados usadas:

| Dado | Arquivo |
|---|---|
| Alimentos candidatos e preços (US$/kg MS) | `GO_campo_dourado/inputs_resolvidos/feed.json`, ração `lac_cow` |
| Composição dos alimentos (base matéria seca) | `GO_campo_dourado/inputs_resolvidos/NASEM_Comp.csv` (`rufas_id` = `feed_type`) |
| Exigências do lote (lado direito das restrições) | CSV de saída do RuFaS, colunas `avg_rqmts_for_LAC_COW_PEN_3.*` |
| Limites de formulação (fibra, gordura, proteína, CMS) | `RUFAS/biophysical/animal/animal_module_constants.py` |
| Restrições usadas pelo otimizador do RuFaS | `RUFAS/biophysical/animal/ration/ration_optimizer.py`, `set_constraints` |

Todos os números de saída do RuFaS citados aqui vêm da **rodada oficial de 5 anos** da fazenda
(`GO_campo_dourado/_pipeline_task_manager_metadata.json`, executada em 25/10/2026; 41 formulações
do lote de lactação). Como rodar e ler o CSV: seção 11.

---

## Bloco 1 — Estágio do animal

O que define a vaca "média" do lote. No quadro: **BW 750 kg, 2ª lactação, DEL 130 dias, 40 kg de leite/dia**.
No RuFaS, o lote de lactação (`LAC_COW_PEN_3`) tem um animal médio recalculado a cada 45 dias
(`ration.formulation_interval`); o CSV reporta `avg_BW` e as exigências desse animal médio.

| Parâmetro | Símbolo | Quadro do professor | Campo Dourado (RuFaS, lote de lactação) |
|---|---|---|---|
| Peso vivo | BW | 750 kg | 635 kg (média do lote, `avg_BW`) |
| Ordem de lactação | — | 2ª | mistura (36 % 1ª, 26 % 2ª, 18 % 3ª, 20 % 4ª+) |
| Dias em lactação | DEL | 130 | varia por vaca; o lote mistura estágios |
| Produção de leite | L | 40 kg/dia | 40,0 kg/vaca/dia (média dos 5 anos) |
| Gordura / proteína do leite | G / Pr | — | 4,0 % / 3,2 % (`animal.json`) |
| Ordenhas por dia | — | — | 3 |

A fazenda tem 4 lotes (o círculo com 15, 20, 30, 35 no quadro): bezerras (`CALF_PEN_0`),
novilhas (`GROWING_PEN_1`), pré-parto (`CLOSE_UP_PEN_2`) e lactação (`LAC_COW_PEN_3`).
**Cada lote é um problema de PL separado**, com seus próprios alimentos permitidos e exigências.
O trabalho foca no lote de lactação, que concentra o custo de alimentação.

---

## Bloco 2 — Requisitos nutricionais (os b_i)

O padrão nutricional (`config.nutrient_standard` = NASEM) transforma o estágio do animal em
exigências diárias. São os **lados direitos** das restrições. Há duas formas de obtê-los:

### 2.1 Direto do RuFaS (forma recomendada)

Rodar a simulação e ler no CSV as colunas `avg_rqmts_for_LAC_COW_PEN_3.*` (uma linha por formulação).
Médias das 41 formulações da rodada oficial:

| Exigência | Coluna no CSV | Média | Faixa nos 5 anos | Unidade |
|---|---|---|---|---|
| EL de mantença | `NEmaint_requirement` | 12,59 | 12,4 a 12,8 | Mcal/dia |
| EL de atividade | `NEa_requirement` | 0,00 | — | Mcal/dia |
| EL de crescimento | `NEg_requirement` | 0,67 | 0,6 a 0,7 | Mcal/dia |
| EL de gestação | `NEpreg_requirement` | 2,20 | 1,4 a 3,7 | Mcal/dia |
| EL de lactação | `NEl_requirement` | 30,01 | 28,4 a 32,8 | Mcal/dia |
| **Energia líquida total** | soma das cinco | **45,5** | | Mcal/dia |
| Proteína metabolizável | `MP_requirement` | **2 575** | 2 404 a 2 855 | g/dia |
| Cálcio absorvível | `Ca_requirement` | **65,2** | 62,1 a 69,5 | g/dia |
| Fósforo absorvível | `P_req` (RuFaS usa o maior entre `P_req` e `P_req_process`) | **61,7** | 59,5 a 65,7 | g/dia |
| Consumo de MS estimado | `DMIest_requirement` | **24,4** | 23,9 a 25,3 | kg/dia |
| Peso vivo médio | `avg_BW` | 635 | 621 a 646 | kg |

A faixa mostra que as exigências mudam a cada formulação (clima, estágio médio do lote). O modelo
usa a média; a análise de sensibilidade nos b_i cobre a variação.

Duas execuções da mesma simulação (mesma seed 44) deram médias cerca de 2 % diferentes. A causa não
foi identificada (o chute inicial do otimizador, `_build_initial_value`, usa `random`, que o RuFaS semeia
em `task_manager.py`); a dieta escolhida muda, e isso realimenta peso e leite do rebanho.
Mais um motivo para tratar os b_i como valores com faixa, não como constantes exatas.

### 2.2 À mão, pelas equações NASEM que o RuFaS implementa

Útil para a vaca do quadro (750 kg, 40 kg de leite, 4 % gordura, 3,2 % proteína). Fórmulas de
`RUFAS/biophysical/animal/nutrients/nasem_requirements_calculator.py` (termos de gestação omitidos):

| Exigência | Equação (NASEM 2021, como no RuFaS) | Vaca do quadro |
|---|---|---|
| EL mantença | 0,10 × BW^0,75 | 0,10 × 750^0,75 = 14,3 Mcal |
| EL lactação | L × (0,0929 G + 0,0588 Pr + 0,0395 Lac), Lac ≈ 4,85 % | 40 × 0,751 = 30,1 Mcal |
| CMS previsto (2ª lactação ou mais) | [9,4 + 0,305 ELl + 0,022 BW − 2,56 ECC] × [1 − 0,348 e^(−0,053 DEL)], ECC = 3 | ≈ 27,4 kg |
| Ca absorvível | 0,90 CMS + L (0,295 + 0,239 Pr) | 24,6 + 42,4 = 67 g |
| P absorvível | CMS + 0,0006 BW + L (0,49 + 0,13 Pr) | 27,9 + 36,2 = 64 g |
| Proteína metabolizável | [0,20 BW^0,6 × 0,85 + 0,73 (11,62 + 0,134 FDN%) CMS + 10 Pr L] / 0,69 + 0,331 BW | ≈ 2 570 g |

Os valores batem em ordem de grandeza com o que o RuFaS reporta para o lote (a vaca do lote é
mais leve, 635 kg, por isso mantença menor; a proteína do lote sai um pouco maior porque o RuFaS
inclui gestação e crescimento).

### 2.3 Limites que não dependem do animal (constantes do RuFaS)

| Limite | Valor | Constante |
|---|---|---|
| FDN total da dieta | entre 25 % e 45 % da MS | `MINIMUM_RATION_NDF`, `MAXIMUM_RATION_NDF` |
| FDN vinda de forragem | ≥ 15 % da MS | `MINIMUM_RATION_FORAGE_NDF` |
| Gordura (extrato etéreo) | ≤ 7 % da MS | `MAXIMUM_RATION_FAT` |
| Proteína metabolizável | ≤ 1,5 × exigência | `PROTEIN_UPPER_LIMIT_FACTOR` |
| Consumo de MS | entre 70 % e 130 % do CMS estimado (otimização livre) | `DMI_CONSTRAINT_FRACTION = 0,30` |

Sobre o **buffer de fósforo** (`ration.phosphorus_requirement_buffer = 50`, citado no PDF do professor
como "exigir 150 % do P"): nesta versão do RuFaS (1.0.5) o parâmetro só existe no esquema de validação
da entrada (`RUFAS/input/metadata/properties/default.json`); nenhum trecho do código o lê. O `P_req`
do CSV (61,7 g) é o requisito NASEM puro, sem buffer. No modelo de PL o buffer entra como um
multiplicador do lado direito (seção 5.4), e é um bom candidato para a análise de sensibilidade,
mas não se deve esperar que o RuFaS reaja a ele.

---

## Bloco 3 — Otimização (o modelo de PL)

### 3.1 Variáveis de decisão

x_j = quantidade do alimento j na dieta, em **kg de matéria seca por vaca por dia**, j = 1..n.
Alimentos permitidos para o lote de lactação (ração `lac_cow` do `feed.json`):

| j | Alimento | `feed_type` | Categoria | MS (%) | c_j (US$/kg MS) |
|---|---|---|---|---|---|
| 1 | Silagem de milho | 51 | Forragem | 35,4 | 0,16 |
| 2 | Silagem de leguminosa | 110 | Forragem | 42,9 | 0,18 |
| 3 | Milho grão moído | 44 | Concentrado | 86,9 | 0,19 |
| 4 | Farelo de soja 48 % | 167 | Concentrado | 89,3 | 0,41 |
| 5 | Caroço de algodão | 58 | Concentrado | 91,4 | 0,36 |
| 6 | Núcleo mineral | 301 | Mineral | 99,5 | 0,80 |

Por que matéria seca e não "kg no cocho": as exigências NASEM e a tabela de composição são em MS.
Para converter depois: kg de matéria natural = x_j / (MS_j / 100). Ex.: 10 kg MS de silagem de
milho = 10 / 0,354 = 28,3 kg no cocho. Se o preço de mercado vier em R$/ton de matéria natural,
c_j (por kg MS) = preço_MN / 1000 / (MS_j / 100).

### 3.2 Função objetivo

Minimizar o custo diário da dieta de uma vaca:

```
min Z = Σ_j c_j x_j = 0,16 x1 + 0,18 x2 + 0,19 x3 + 0,41 x4 + 0,36 x5 + 0,80 x6     [US$/vaca/dia]
```

É a mesma função que o RuFaS minimiza (`RationOptimizer.objective`: soma de quantidade × preço).
O custo do lote é Z × número de vacas; o da fazenda, a soma dos lotes.

### 3.3 Forma geral das restrições

Cada nutriente i gera uma restrição linear "o que a dieta entrega" versus "o que a vaca precisa":

```
Σ_j a_ij x_j  ≥  b_i     (mínimos: energia, proteína, Ca, P)
Σ_j a_ij x_j  ≤  b_i     (máximos: proteína, CMS)
```

onde a_ij = quantidade do nutriente i em 1 kg de MS do alimento j, e b_i = exigência do Bloco 2.

Restrições dadas em **porcentagem da dieta** (FDN, FDN de forragem, gordura) são razões
Σ a_ij x_j / Σ x_j, que não são lineares. Multiplicando pelo CMS (Σ x_j > 0) ficam lineares:

```
Σ_j a_ij x_j / Σ_j x_j ≥ k    ⇔    Σ_j (a_ij − k) x_j ≥ 0
Σ_j a_ij x_j / Σ_j x_j ≤ k    ⇔    Σ_j (k − a_ij) x_j ≥ 0
```

---

## 4. Coeficientes técnicos a_ij (por kg de MS)

Da tabela NASEM. Energia e proteína metabolizável foram **linearizadas** (detalhes na seção 6).

| Nutriente i (por kg MS) | x1 Sil. milho | x2 Sil. legum. | x3 Milho | x4 F. soja | x5 Car. algodão | x6 Mineral |
|---|---|---|---|---|---|---|
| Energia líquida, Mcal | 1,75 | 1,54 | 2,06 | 2,36 | 1,87 | 0 |
| Proteína metabolizável, g | 77,8 | 92,1 | 96,8 | 236,5 | 102,6 | 4,7 |
| Ca absorvível, g | 0,71 | 3,76 | 0,21 | 2,40 | 1,04 | 146,6 |
| P absorvível, g | 1,48 | 2,25 | 2,16 | 5,17 | 4,35 | 18,56 |
| FDN, % | 40,93 | 43,20 | 9,76 | 11,06 | 50,58 | 0,14 |
| FDN de forragem, % | 40,93 | 43,20 | 0 | 0 | 0 | 0 |
| Extrato etéreo (gordura), % | 2,92 | 2,87 | 3,85 | 1,81 | 18,62 | 0 |

Como foram obtidos (colunas de `NASEM_Comp.csv`):

- **Ca e P absorvíveis** = teor (%) × 10 × digestibilidade que o RuFaS aplica por categoria
  (`calculate_calcium_supply` / `calculate_phosphorus_supply`): Ca forragem 0,30, concentrado 0,60,
  mineral 0,95; P forragem 0,64, concentrado 0,70, mineral 0,80.
- **Energia líquida** = k × `DE_Base`, com k = 0,596. No RuFaS, EL = 0,66 × EM (`EFF_OF_ME_USE`) e
  EM = ED da dieta − energia do metano − energia do N urinário. Na dieta do caso base, EM/ED = 0,90,
  logo k ≈ 0,66 × 0,90 ≈ 0,60; o valor 0,596 foi ajustado para reproduzir a EL que o RuFaS reportou
  para a dieta praticada (45,6 Mcal; com 0,60 exato daria 45,9). Se usasse a relação EM/ED "de livro" (0,82), k = 0,54 e o modelo
  subestimaria a energia em 9 %, forçando dietas mais caras do que o necessário.
- **Proteína metabolizável** = 0,64 × 0,13 × NDT(g/kg) + PB(g/kg) × PNDR × dPNDR + 4,72 (endógena),
  com NDT = `DE_Base` / 0,04409, PNDR = coluna `RUP`, dPNDR = coluna `dRUP`. Reproduz a conta do RuFaS
  supondo que a proteína microbiana seja limitada pelo NDT, o que ocorre em dietas típicas de lactação
  (erro de 2 % na dieta do caso base, seção 9).
- **FDN, FDN de forragem, EE**: colunas `NDF` e `EE`; FDN de forragem só conta para `feed_type = Forage`.

---

## 5. Restrições do modelo (instância Campo Dourado)

Notação: CMS = x1 + x2 + x3 + x4 + x5 + x6.

### 5.1 Energia líquida (≥ exigência total)
```
1,75 x1 + 1,54 x2 + 2,06 x3 + 2,36 x4 + 1,87 x5 + 0 x6  ≥  45,5
```

### 5.2 Proteína metabolizável (faixa)
```
77,8 x1 + 92,1 x2 + 96,8 x3 + 236,5 x4 + 102,6 x5 + 4,7 x6  ≥  2 575        (mínimo)
77,8 x1 + 92,1 x2 + 96,8 x3 + 236,5 x4 + 102,6 x5 + 4,7 x6  ≤  3 863        (1,5 × exigência)
```

### 5.3 Cálcio absorvível (≥)
```
0,71 x1 + 3,76 x2 + 0,21 x3 + 2,40 x4 + 1,04 x5 + 146,6 x6  ≥  65,2
```

### 5.4 Fósforo absorvível (≥)
```
1,48 x1 + 2,25 x2 + 2,16 x3 + 5,17 x4 + 4,35 x5 + 18,56 x6  ≥  61,7 × (1 + buffer/100)
```
Caso base: buffer = 0 → 61,7 g. Cenário "buffer 50 %": 92,6 g.

### 5.5 FDN total entre 25 % e 45 % da MS
```
mínimo:  15,93 x1 + 18,20 x2 − 15,24 x3 − 13,94 x4 + 25,58 x5 − 24,86 x6  ≥  0
máximo:   4,07 x1 +  1,80 x2 + 35,24 x3 + 33,94 x4 −  5,58 x5 + 44,86 x6  ≥  0
```

### 5.6 FDN de forragem ≥ 15 % da MS
```
25,93 x1 + 28,20 x2 − 15 x3 − 15 x4 − 15 x5 − 15 x6  ≥  0
```

### 5.7 Gordura (extrato etéreo) ≤ 7 % da MS
```
4,08 x1 + 4,13 x2 + 3,15 x3 + 5,19 x4 − 11,62 x5 + 7,00 x6  ≥  0
```

### 5.8 Consumo de matéria seca (faixa em torno do CMS estimado)
```
x1 + x2 + x3 + x4 + x5 + x6  ≥  0,70 × 24,4 = 17,0
x1 + x2 + x3 + x4 + x5 + x6  ≤  1,30 × 24,4 = 31,7
```

Por que faixa e não igualdade: é o que o RuFaS faz na otimização livre (±30 %), e o CMS estimado é
uma previsão, não uma medida. Mas atenção ao efeito no Simplex: como todo c_j > 0 e a energia é a
restrição que mais pesa, o modelo escolhe livremente o CMS que fecha a energia ao menor custo, e pode
sair acima do estimado (a vaca "come mais de comida barata"). Duas variantes úteis para comparar:

- CMS fixo: trocar a faixa por Σ x_j = 24,4 (uma igualdade, que no Simplex entra com variável artificial).
  Com os coeficientes acima a instância continua viável.
- Faixa estreita: ±10 % (no RuFaS com `user_input = true` a tolerância de 10 % é por alimento, em torno
  de ρ_j × 1,1 × CMS, o que na soma dá um CMS entre 0,99 e 1,21 do estimado).

### 5.9 Limites por alimento e não-negatividade
```
L_j ≤ x_j ≤ U_j,   x_j ≥ 0,   j = 1..6
```
No RuFaS, L_j e U_j vêm das colunas `lower_limit` (0) e `limit` (100 kg) da tabela de composição:
na prática, sem limite. Isso importa: uma solução básica do Simplex tem no máximo tantas variáveis
positivas quanto restrições ativas, então sem limites de inclusão a dieta ótima tende a se concentrar
em dois ou três alimentos e zerar o resto (mineral a 0 kg, por exemplo). Para uma dieta praticável
vale acrescentar limites agronômicos (ex.: caroço de algodão ≤ 3 kg por causa do gossipol; mineral
≥ 0,2 kg; silagem de milho ≥ 20 % da MS). Cada limite é uma restrição a mais no Simplex (ou um
bound, se o solver aceitar).

---

## 6. O que é exato e o que é aproximação em relação ao RuFaS

O RuFaS resolve o mesmo problema com um método **não linear** (SLSQP, `scipy.optimize.minimize`),
porque duas de suas restrições não são lineares em x:

| Item | RuFaS (não linear) | Este modelo (PL) |
|---|---|---|
| Função objetivo | Σ c_j x_j | idêntica |
| Ca, P, CMS | lineares | idênticas |
| FDN, FDN forragem, gordura | razões Σ a x / Σ x | reescritas como Σ (a − k) x ≥ 0 (equivalência exata para Σ x > 0) |
| Energia | ED NASEM de cada alimento depende do CMS total, do peso vivo e do amido total da dieta (`calculate_NASEM_dNDF`, `calculate_NASEM_dstarch`); depois desconta metano e N urinário do lote | coeficiente fixo por alimento, k × DE_Base, calibrado na dieta do caso base |
| Proteína metabolizável | min(0,13 NDT; 0,85 PDR) × 0,64 + PNDR + endógena, com taxas de passagem que dependem do CMS/BW e do % de concentrado | assume limitação por NDT e PNDR fixa da tabela |
| Limite superior de PM | 1,5 × exigência | igual |

Em PO isso se chama **linearização em torno de um ponto de operação**: os coeficientes de energia e
proteína valem para dietas parecidas com a da fazenda (CMS ≈ 25 kg, ~50 % forragem). Se a solução
ótima do Simplex ficar longe desse ponto, recalculam-se os coeficientes e resolve-se de novo
(uma iteração de "PL sequencial"). A validação final é rodar o RuFaS com a dieta do Simplex
(`ration.user_input = true`) e conferir no CSV `ration_nutrient_amount_for_LAC_COW_PEN_3.*`.

Vantagens do Simplex sobre o SLSQP, para o relatório: ótimo **global** garantido (o SLSQP é local e
parte de um chute inicial, aleatório na primeira formulação, `_build_initial_value`), e **preços-sombra, custos reduzidos e faixas de
sensibilidade** saem de graça do quadro final.

---

## 7. Modelo completo e forma padrão para o Simplex

### 7.1 Forma compacta

```
min  Z = Σ_j c_j x_j

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
     x_j ≥ 0
```

6 variáveis, 11 restrições funcionais.

### 7.2 Matriz do problema (para o tableau ou para o Solver)

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

### 7.3 Forma padrão

- Restrições "≤" (R3, R11): soma-se uma variável de **folga** s_i ≥ 0.
- Restrições "≥" (R1, R2, R4..R10): subtrai-se uma variável de **excesso** e_i ≥ 0 e soma-se uma
  **artificial** a_i ≥ 0, porque a origem x = 0 não é viável (energia, proteína, Ca, P e CMS mínimo falham).
- Resolve-se por **Duas Fases** (Fase 1 minimiza Σ a_i; Fase 2 minimiza Z) ou **Big-M**.

Tamanho do tableau: 11 linhas × (6 x + 2 folgas + 9 excessos + 9 artificiais) = 11 × 26 colunas.
Como as linhas R6..R9 têm b = 0, a solução pode ser degenerada; a regra de Bland evita ciclagem.

Dica de escala: as linhas de proteína têm coeficientes na casa das centenas e a de energia na casa
das unidades. Para fazer à mão, vale dividir R2 e R3 por 100 (PM em "centenas de g"), o que não muda
a solução, só a unidade do preço-sombra.

### 7.4 Leitura do quadro ótimo

- Variáveis básicas → x_j ótimos (kg MS). Dividir por CMS → `ration_percentage` para o `feed.json`.
- Linha Z → custo mínimo (US$/vaca/dia).
- Coeficientes das folgas/excessos na linha Z → **preços-sombra** y_i: quanto o custo sobe se a
  exigência b_i aumentar 1 unidade (ex.: US$ por Mcal de energia, por g de P).
- Restrições com excesso/folga = 0 são **ativas** (gargalos); as demais são folgadas.
- Custos reduzidos das não básicas: quanto o preço de um alimento ausente precisa cair para entrar na dieta.

---

## 8. Extensões lineares (ficam para depois, mas cabem no Simplex)

- **Vários lotes com estoque comum**: x_pj por lote p, n_p vacas por lote;
  min Σ_p n_p Σ_j c_j x_pj, com as restrições nutricionais de cada lote e
  Σ_p n_p x_pj ≤ Disp_j (kg/dia disponíveis de silagem, `*_purchase_allowance`). Acopla os lotes
  e continua linear.
- **Horizonte de tempo**: um bloco de variáveis por intervalo de formulação (8 por ano), com estoque
  de silagem transportado entre períodos. É um PL maior, não um problema novo.
- **Efeito econômico**: margem = preço do leite × leite − Z. O RuFaS não tem preço de leite
  (`economy.csv` só traz diesel, energia, gás e água); entra por fora. Na PL, a produção é dado
  (exigência), não variável; para otimizar leite seria preciso tornar L variável, e a exigência de
  energia passa a depender de L (ainda linear: ELl = 0,751 L).

---

## 9. Checagem do modelo contra o caso base do RuFaS

Dieta que o RuFaS praticou no caso base (média das 41 formulações, `ration_per_animal_for_LAC_COW_PEN_3.*`):
silagem de milho 7,87, silagem de leguminosa 4,34, milho 5,80, farelo de soja 3,88, caroço de algodão 2,17,
mineral 0,72 kg MS (CMS 24,8 kg; em % de MS: 32 / 18 / 23 / 16 / 9 / 3). Avaliando essa dieta com os
coeficientes da seção 4 e comparando com o que o próprio RuFaS reportou em `ration_nutrient_amount_for_LAC_COW_PEN_3.*`:

| Nutriente | Modelo PL | RuFaS reportou | Exigência / limite | Situação |
|---|---|---|---|---|
| Custo, US$/vaca/dia | 6,09 | (não reporta; mesma conta) | — | |
| Energia líquida, Mcal | 45,6 | 45,6 | ≥ 45,5 | **ativa** (gargalo) |
| Proteína metabolizável, g | 2 716 | 2 669 | 2 575 a 3 863 | folga pequena |
| Ca absorvível, g | 141 | 141 | ≥ 65,2 | grande folga (mineral) |
| P absorvível, g | 77 | 77 | ≥ 61,7 | folga |
| FDN, % MS | 29,0 | — | 25 a 45 | folga |
| FDN de forragem, % MS | 20,6 | 20,6 | ≥ 15 | folga |
| Gordura, % MS | 4,24 | 4,24 | ≤ 7 | folga |

Conclusões:

- Os coeficientes linearizados reproduzem o RuFaS com erro ≤ 2 % em todos os nutrientes; Ca, P, FDN
  de forragem e gordura batem exatamente (as colunas de Ca e P do CSV já são a fração absorvível).
- No caso base só a **energia** está ativa. É ela que define o custo; o preço-sombra da energia
  (US$/Mcal) será o número mais importante do relatório de sensibilidade.
- O mineral a 3 % entrega o dobro do cálcio necessário a US$ 0,80/kg: há custo a cortar, e é o
  primeiro lugar onde o Simplex vai mexer.
- A instância é **viável** (verifiquei resolvendo-a com e sem CMS fixo). O valor ótimo serve de gabarito
  para conferir o Simplex feito à mão e, de propósito, não está neste documento.

---

## 10. Próximos passos do projeto (depois da modelagem)

1. Resolver a instância da seção 7 por Simplex (à mão para o entendimento; Solver do LibreOffice/Excel
   em modo "Programação linear" para conferir) e pedir o relatório de sensibilidade.
2. Converter x_j em `ration_percentage`, colocar no `feed.json` com `ration.user_input = true`
   e rodar o RuFaS (seção 11) para validar nutrientes e custo.
3. Efeito econômico: custo/vaca/dia × vacas em lactação × 365, antes e depois; margem com preço de leite externo.
4. Sensibilidade: os 7 cenários do PDF do professor (preço do milho, preço das silagens, buffer de P,
   NASEM × NRC, intervalo de formulação, CMS/leite ±10 %, limites de compra), comparando os pontos de
   quebra observados com as faixas que o Simplex previu.
5. Só então: o módulo de otimização integrado ao RuFaS.

---

## 11. Como reproduzir os números deste documento

Tudo parte da raiz do RuFaS e do pacote do professor (`material_complementar/Fazenda Campo Dourado (1).zip`,
extraído em `GO_campo_dourado/`). Nenhum arquivo-fonte do RuFaS é alterado.

**Ambiente** (uma vez). O RuFaS exige Python 3.12 ou 3.13; o `python3` padrão do Mac é 3.14:

```bash
/opt/homebrew/bin/python3.13 -m venv venv
venv/bin/pip install .
```

**Filtro de saída**. O RuFaS grava no CSV só as variáveis que casam com os padrões dos arquivos
`output/output_filters/csv_*.txt`. O arquivo `csv_dieta.txt` (6 linhas, criado para este trabalho)
seleciona ração por animal, exigências médias, nutrientes entregues, ração do usuário e leite.
Sem ele, o CSV sai com milhares de colunas.

**Simulação oficial** (5 anos, fazenda completa, cerca de 2 minutos):

```bash
venv/bin/python main.py -g -p GO_campo_dourado/_pipeline_task_manager_metadata.json
```

Deve terminar com `Finished task: 1/1 (GO_campo_dourado) with 0 error(s)`. O resultado fica em
`output/CSVs/GO_campo_dourado_saved_variables_csv_dieta.txt_<data>.csv`. Abra numa planilha.

**Leitura do CSV.** Cada coluna é `Classe.função.variável.subchave (unidade)`. As linhas não são
dias: cada variável ocupa só as primeiras linhas de que precisa (41 para as variáveis de ração, uma
por formulação; 1 826 para leite, uma por dia). Colunas do lote de lactação (`LAC_COW_PEN_3`):

| Grupo de colunas | O que é | Uso no modelo |
|---|---|---|
| `avg_rqmts_for_LAC_COW_PEN_3.*` | exigências do animal médio do lote | lados direitos b_i (seção 2.1) |
| `ration_per_animal_for_LAC_COW_PEN_3.<feed_type> (kg)` | kg de MS de cada alimento por vaca/dia | a dieta praticada (seção 9); custo = Σ kg × preço do `feed.json` |
| `ration_nutrient_amount_for_LAC_COW_PEN_3.*` | o que a dieta entregou (EL, PM, Ca, P, FDN de forragem, gordura) | conferência dos coeficientes a_ij (seção 9) |
| `user_defined_ration.<feed_type> (percent)` | percentagens escritas no `feed.json` | |
| `daily_milk_production (kg/day)`, `milking_cow_num (animals)` | leite do rebanho e vacas em ordenha, por dia | leite por vaca = razão entre os dois |

Atenção: as colunas `calcium` e `phosphorus` de `ration_nutrient_amount` já são a fração absorvível
(digestibilidade aplicada), em kg; multiplique por 1 000 para comparar com `Ca_requirement` e `P_req`, em g.
