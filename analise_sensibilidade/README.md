# Análise de sensibilidade

Painel Streamlit sobre o tableau ótimo que o Simplex do RuFaS grava em `output/simplex/GO_campo_dourado/`.
Sliders de preço e de exigência, re-otimização do PL, tableau e faixas de validade de c_j e b_i.

Rodar, da raiz do RuFaS (precisa ter rodado `main.py -g -r simplex -p GO_campo_dourado/...` antes):

```
venv/bin/pip install -r analise_sensibilidade/requirements.txt
venv/bin/streamlit run analise_sensibilidade/app.py
```

No Windows, `venv\Scripts\...`. Só o venv, sem Docker. O `.streamlit/config.toml` na raiz desliga o pedido de e-mail.

Arquivos: `app.py` (painel), `sensibilidade.py` (leitura dos CSVs, faixas de c_j e b_i, re-solução com o
solver do RuFaS), `RESULTADOS.md` (números para o relatório).

Trabalhos futuros: alimento novo (custo reduzido de coluna nova), limite de kg por alimento, mudar
composição a_ij, levar a dieta simulada de volta ao RuFaS.
