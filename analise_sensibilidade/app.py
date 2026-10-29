"""
Painel de análise de sensibilidade do Simplex do RuFaS (Faz. Campo Dourado).

Lê o tableau ótimo gravado pelo RuFaS, deixa o usuário mexer nos preços dos alimentos e nas exigências
da vaca com sliders, re-otimiza o PL e mostra a dieta nova, o custo, o tableau e as faixas de validade.

Executar da raiz do RuFaS:  venv/bin/streamlit run analise_sensibilidade/app.py
"""

from __future__ import annotations

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

import sensibilidade as sens

AZUL = "#2a78d6"  # dieta atual
LARANJA = "#eb6834"  # dieta simulada

st.set_page_config(page_title="Análise de sensibilidade", layout="wide")

# Exigências que viram slider. Fibra, gordura e consumo mínimo são constantes nutricionais e ficam fixas.
EXIGENCIAS = {
    "NASEM_net_energy_constraint": ("Energia líquida", "Mcal/dia"),
    "protein_constraint_lower": ("Proteína metabolizável", "g/dia"),
    "calcium_constraint": ("Cálcio", "g/dia"),
    "phosphorus_constraint": ("Fósforo", "g/dia"),
    "DMI_constraint_upper": ("Consumo máximo de MS", "kg/dia"),
}


def br(valor: float, casas: int = 2) -> str:
    """Número com vírgula decimal."""
    if valor is None or not np.isfinite(valor):
        return "-"
    return f"{valor:,.{casas}f}".replace(",", "X").replace(".", ",").replace("X", ".")


def rotulo_formulacao(f: sens.Formulacao) -> str:
    ano, dia = divmod(f.day, 365)
    return f"Dia {f.day} da simulação (ano {ano + 1}, dia {dia})"


def restaurar_padroes(padroes: dict[str, float]) -> None:
    st.session_state.update(padroes)


def grafico_dietas(atual: pd.DataFrame, nova: pd.DataFrame) -> alt.Chart:
    dados = pd.concat([atual.assign(dieta="Dieta atual"), nova.assign(dieta="Dieta simulada")]).reset_index()
    return (
        alt.Chart(dados)
        .mark_bar(cornerRadiusEnd=4)
        .encode(
            x=alt.X("kg no cocho:Q", title="kg por vaca por dia (no cocho)"),
            y=alt.Y("alimento:N", title=None, sort=None),
            yOffset=alt.YOffset("dieta:N"),
            color=alt.Color(
                "dieta:N",
                title=None,
                scale=alt.Scale(domain=["Dieta atual", "Dieta simulada"], range=[AZUL, LARANJA]),
                legend=alt.Legend(orient="top"),
            ),
            tooltip=[
                alt.Tooltip("alimento:N", title="Alimento"),
                alt.Tooltip("dieta:N", title="Dieta"),
                alt.Tooltip("kg no cocho:Q", title="kg no cocho", format=".2f"),
                alt.Tooltip("kg MS/vaca/dia:Q", title="kg MS", format=".2f"),
                alt.Tooltip("% da MS:Q", title="% da MS", format=".1f"),
            ],
        )
        .properties(height=280, width="container")
    )


def tabela_dieta(modelo: sens.Modelo, x: np.ndarray, dm: dict[str, float]) -> pd.DataFrame:
    tabela = sens.dieta(modelo, x)
    tabela["kg no cocho"] = [v / (dm[f] / 100.0) for v, f in zip(tabela["kg MS/vaca/dia"], modelo.feeds)]
    return tabela


def nome_legivel(nome: str) -> str:
    """Nomes das linhas e colunas do tableau."""
    if nome in ("Z", "RHS"):
        return nome
    if nome in sens.FEED_NAMES:
        return sens.FEED_NAMES[nome]
    for prefixo, rotulo in (("folga_max_", "folga teto "), ("excesso_", "excesso "), ("folga_", "folga ")):
        if nome.startswith(prefixo):
            resto = nome[len(prefixo):]
            if prefixo == "folga_max_":
                return rotulo + sens.FEED_NAMES.get(resto, resto)
            return rotulo + sens.constraint_label(resto).split(" (")[0]
    return nome


def mostrar_tableau(df: pd.DataFrame, base: list[str], z: float) -> None:
    legivel = df.copy()
    legivel.index = [nome_legivel(i) for i in legivel.index]
    legivel.columns = [nome_legivel(c) for c in legivel.columns]
    estilo = legivel.style.format("{:.4f}").set_properties(subset=pd.IndexSlice[["Z"], :], **{"font-weight": "bold"})
    st.dataframe(estilo, width="stretch", height=min(40 + 35 * len(legivel), 700))
    st.write(f"**Z = US$ {br(z, 4)} /vaca/dia** - base ótima: " + ", ".join(nome_legivel(b) for b in base))


# ------------------------------------------------------------------------------------------------
# Barra lateral: escolha da formulação e flexibilidade da receita
# ------------------------------------------------------------------------------------------------

st.sidebar.title("Análise de sensibilidade")
st.sidebar.caption("Fazenda Campo Dourado - lote de vacas em lactação")

fazendas = sens.listar_fazendas()
if not fazendas:
    st.error(
        f"Não há resultados do Simplex em `{sens.SIMPLEX_DIR}`. Rode antes, da raiz do RuFaS:\n\n"
        "```bash\nvenv/bin/python main.py -g -r simplex -p GO_campo_dourado/_pipeline_task_manager_metadata.json\n```"
    )
    st.stop()

fazenda = fazendas[0] if len(fazendas) == 1 else st.sidebar.selectbox("Simulação", fazendas, format_func=lambda p: p.name)
formulacoes = sens.listar_formulacoes(fazenda)
formulacao = st.sidebar.selectbox("Momento da simulação", formulacoes, index=len(formulacoes) - 1, format_func=rotulo_formulacao)
modelo = sens.ler_modelo(formulacao.modelo)
tableau = sens.ler_tableau(formulacao.tableau)
dm = sens.materia_seca(modelo.feeds)
custos = sens.faixa_custos(tableau, modelo)
rhs = sens.faixa_rhs(tableau, modelo)
funcionais = list(modelo.funcionais.index)
c0 = modelo.variaveis.loc["custo_c", modelo.feeds].to_numpy(dtype=float)
b0 = modelo.funcionais["b"].to_numpy(dtype=float)
x_atual = modelo.variaveis.loc["x_otimo", modelo.feeds].to_numpy(dtype=float)
L0 = modelo.variaveis.loc["limite_inferior", modelo.feeds].to_numpy(dtype=float)
U0 = modelo.variaveis.loc["limite_superior", modelo.feeds].to_numpy(dtype=float)

# Valores iniciais dos sliders (os da formulação escolhida).
padroes = {"flex": 10} | {f"preco_{f}": float(c0[j]) for j, f in enumerate(modelo.feeds)}
for nome in EXIGENCIAS:
    if nome in funcionais:
        padroes[f"exig_{nome}"] = float(rhs.loc[rhs["nome"] == nome, "b_exib"].iloc[0])
if st.session_state.get("formulacao_atual") != formulacao.label:
    st.session_state["formulacao_atual"] = formulacao.label
    restaurar_padroes(padroes)
flex = st.sidebar.slider(
    "Quanto a dieta pode se afastar da receita",
    min_value=10,
    max_value=100,
    step=10,
    format="±%d %%",
    key="flex",
    help="10 % é o que o RuFaS usa: cada alimento fica entre 90 % e 110 % da receita da fazenda.",
)

st.sidebar.button("Voltar aos valores atuais", width="stretch", on_click=restaurar_padroes, args=(padroes,))

# ------------------------------------------------------------------------------------------------
# Página
# ------------------------------------------------------------------------------------------------

st.title("Análise de sensibilidade da dieta de menor custo")
st.write(
    "Ajuste os preços dos alimentos e as exigências da vaca. A cada mudança o Simplex re-otimiza a dieta "
    "e o painel mostra o que muda em relação à dieta atual, o tableau e as faixas de validade."
)

# --- 1. Preços -----------------------------------------------------------------------------------
st.subheader("1. Preços dos alimentos (US$/kg MS)")
colunas = st.columns(3)
novos_precos = c0.copy()
for j, feed in enumerate(modelo.feeds):
    with colunas[j % 3]:
        novos_precos[j] = st.slider(
            sens.FEED_NAMES.get(feed, f"Alimento {feed}"),
            min_value=0.01,
            max_value=float(round(max(c0[j] * 2.5, 0.10), 2)),
            step=0.01,
            format="US$ %.2f",
            key=f"preco_{feed}",
            help=f"Preço atual: US$ {c0[j]:.2f}/kg MS.",
        )

# --- 2. Exigências -------------------------------------------------------------------------------
st.subheader("2. Exigências da vaca por dia")
st.caption("Valores calculados pelo RuFaS (NASEM 2021) para a vaca média do lote neste dia.")
novos_b = b0.copy()
exigencias = [n for n in EXIGENCIAS if n in funcionais]
colunas = st.columns(len(exigencias))
for coluna, nome in zip(colunas, exigencias):
    rotulo, unidade = EXIGENCIAS[nome]
    linha = rhs[rhs["nome"] == nome].iloc[0]
    fator, valor_atual = int(linha["fator"]), float(linha["b_exib"])
    with coluna:
        novo = st.slider(
            rotulo,
            min_value=float(round(valor_atual * 0.7, 1)),
            max_value=float(round(valor_atual * 1.3, 1)),
            step=float(max(round(valor_atual / 200, 2), 0.1)),
            format=f"%.1f {unidade}",
            key=f"exig_{nome}",
            help=f"Valor atual: {br(valor_atual, 1)} {unidade}.",
        )
    novos_b[funcionais.index(nome)] = fator * novo
    if nome == "protein_constraint_lower" and "protein_constraint_upper" in funcionais and valor_atual:
        # O teto de proteína acompanha a exigência (1,5 x no RuFaS).
        k = funcionais.index("protein_constraint_upper")
        novos_b[k] = b0[k] * (novo / valor_atual)

# --- 3. Resultado --------------------------------------------------------------------------------
st.subheader("3. Resultado")
receita = (L0 + U0) / 2.0  # centro da faixa ±10 % do RuFaS = receita da fazenda (kg MS)
t = flex / 100.0
lp = modelo.lp(
    c=novos_precos,
    b=novos_b,
    lower=np.maximum(0.0, (1.0 - t) * receita) if flex > 10 else None,
    upper=(1.0 + t) * receita if flex > 10 else None,
)
resultado = sens.resolver(lp)
dieta_atual = tabela_dieta(modelo, x_atual, dm)
custo_manter = float(novos_precos @ x_atual)

if resultado.status != "optimal":
    st.error(
        "Com esses valores não existe dieta que atenda a todas as restrições. "
        "Afrouxe alguma exigência ou aumente a flexibilidade da receita na barra lateral."
    )
else:
    dieta_nova = tabela_dieta(modelo, resultado.x, dm)
    m1, m2, m3 = st.columns(3)
    m1.metric("Custo da dieta atual", f"US$ {br(modelo.Z)} /vaca/dia")
    m2.metric(
        "Custo da dieta simulada",
        f"US$ {br(resultado.objective)} /vaca/dia",
        delta=None if abs(resultado.objective - modelo.Z) < 0.005 else f"{br(resultado.objective - modelo.Z)} US$",
        delta_color="inverse",
    )
    m3.metric(
        "Manter a dieta atual com os preços novos",
        f"US$ {br(custo_manter)} /vaca/dia",
        delta=None if abs(custo_manter - resultado.objective) < 0.005 else f"{br(custo_manter - resultado.objective)} US$ a mais que a simulada",
        delta_color="inverse",
    )

    st.altair_chart(grafico_dietas(dieta_atual, dieta_nova))

    mudancas = [
        (sens.FEED_NAMES.get(f, f), a, b)
        for f, a, b in zip(modelo.feeds, dieta_atual["kg no cocho"], dieta_nova["kg no cocho"])
        if abs(b - a) >= 0.05
    ]
    if mudancas:
        st.info(
            "O que muda na dieta (kg no cocho por vaca por dia):\n"
            + "\n".join(f"- {nome}: {br(a)} -> {br(b)} kg" for nome, a, b in mudancas)
        )
    else:
        st.success("A dieta não muda: com esses valores a receita atual continua sendo a mais barata.")
        presos = [
            sens.FEED_NAMES.get(f, f)
            for j, f in enumerate(modelo.feeds)
            if (novos_precos[j] > c0[j] + 1e-9 and abs(x_atual[j] - lp.lower[j]) < 1e-6)
            or (novos_precos[j] < c0[j] - 1e-9 and abs(x_atual[j] - lp.upper[j]) < 1e-6)
        ]
        if presos:
            st.caption(
                ", ".join(presos) + f" já está no limite da receita (±{flex} %). Para o modelo trocar de alimento, "
                "aumente a flexibilidade na barra lateral."
            )
    if flex >= 50:
        st.warning(
            "Os coeficientes de energia e proteína foram linearizados em torno da receita atual. "
            "Para dietas muito diferentes o resultado é aproximado."
        )

    with st.expander("Dieta em tabela"):
        comparacao = pd.concat(
            [
                dieta_atual[["kg no cocho", "kg MS/vaca/dia", "% da MS"]].add_suffix(" (atual)"),
                dieta_nova[["kg no cocho", "kg MS/vaca/dia", "% da MS"]].add_suffix(" (simulada)"),
            ],
            axis=1,
        )
        st.dataframe(comparacao.style.format("{:.2f}"), width="stretch")

# --- 4. Preços-sombra ----------------------------------------------------------------------------
st.subheader("4. Restrições ativas e preços-sombra da dieta atual")
ativas = rhs[(rhs["tipo"] != "limite_superior") & rhs["ativa"]]
if ativas.empty:
    st.write("Nenhuma exigência nutricional está ativa: a dieta está presa só pelos limites da receita.")
for restricao, linha in ativas.iterrows():
    y = float(linha["preço-sombra_exib"])
    st.write(f"- **{restricao}** está ativa. Preço-sombra: US$ {br(abs(y), 4)} por unidade a mais exigida.")
piso = [sens.FEED_NAMES.get(f, f) for j, f in enumerate(modelo.feeds) if abs(x_atual[j] - L0[j]) < 1e-6]
teto = [sens.FEED_NAMES.get(f, f) for j, f in enumerate(modelo.feeds) if abs(x_atual[j] - U0[j]) < 1e-6]
if piso:
    st.write("- No mínimo permitido pela receita (o modelo usaria menos): " + ", ".join(piso) + ".")
if teto:
    st.write("- No máximo permitido pela receita (o modelo usaria mais): " + ", ".join(teto) + ".")

# --- 5. Tableau ----------------------------------------------------------------------------------
st.divider()
st.header("5. Quadro Simplex (tableau ótimo)")
st.caption(
    "Linhas: variáveis básicas e Z. Colunas: alimentos (variável deslocada x_j - L_j), excesso de cada "
    "restrição >=, folga de cada teto de alimento e o RHS. Na linha Z estão os custos reduzidos e o valor de Z."
)
ajustado = resultado.status == "optimal" and (not np.allclose(novos_precos, c0) or not np.allclose(novos_b, b0) or flex > 10)
if ajustado:
    st.subheader("PL simulado (com os ajustes dos sliders)")
    mostrar_tableau(sens.tableau_para_dataframe(resultado), resultado.basis, resultado.objective)
    st.caption(f"Pivôs do Simplex nesta resolução: {resultado.pivots}.")
    with st.expander("Quadro gravado pelo RuFaS (dieta atual)"):
        mostrar_tableau(tableau, modelo.resumo["base_otima"].split(), modelo.Z)
else:
    st.subheader("Quadro gravado pelo RuFaS (dieta atual)")
    mostrar_tableau(tableau, modelo.resumo["base_otima"].split(), modelo.Z)
    st.caption(
        f"Arquivo `{formulacao.tableau.relative_to(sens.REPO)}` - {modelo.resumo['pivos_simplex']} pivôs - "
        f"{modelo.resumo['iteracoes_linearizacao']} iterações da linearização."
    )

# --- 6. Faixas de validade -----------------------------------------------------------------------
st.divider()
st.header("6. Faixas de validade da base ótima")
st.caption(f"Calculadas do tableau da dieta atual (arquivo `{formulacao.modelo.relative_to(sens.REPO)}`).")

st.subheader("Custos c_j")
st.caption(
    "Dentro de [c_j mínimo, c_j máximo] a base ótima não muda: mesma dieta, Z varia x_j* por unidade de preço. "
    "Alimento não básico: pode subir sem limite e cair até c_j - custo reduzido."
)
st.dataframe(
    custos[["c_j (US$/kg MS)", "x* (kg MS)", "custo reduzido", "situação", "c_j mínimo", "c_j máximo"]]
    .style.format({"c_j (US$/kg MS)": "{:.4f}", "x* (kg MS)": "{:.3f}", "custo reduzido": "{:.4f}"})
    .format(lambda v: br(v, 4), subset=["c_j mínimo", "c_j máximo"]),
    width="stretch",
)

st.subheader("Lados direitos b_i")
st.caption(
    "Dentro de [b mínimo, b máximo] a base ótima não muda e Z varia preço-sombra x delta b. "
    "Máximos mostrados no sentido natural (<=). As linhas 'Limite superior' são os tetos da receita."
)
exib = rhs[["sinal_exib", "b_exib", "valor_exib", "folga", "ativa", "preço-sombra_exib", "b mínimo_exib", "b máximo_exib"]]
exib = exib.rename(columns=lambda c: c.replace("_exib", "")).rename(columns={"valor": "valor em x*"})
st.dataframe(
    exib.style.format({"b": "{:.3f}", "valor em x*": "{:.3f}", "folga": "{:.3f}", "preço-sombra": "{:.4f}"})
    .format(lambda v: br(v, 3), subset=["b mínimo", "b máximo"]),
    width="stretch",
)
