"""
Leitura dos arquivos gravados pelo Simplex do RuFaS (output/simplex/<fazenda>/modelo_*.csv e tableau_*.csv)
e análise de sensibilidade a partir do tableau ótimo: faixas de validade dos custos c_j, faixas dos lados
direitos b_i e nova resolução do PL quando a mudança sai da faixa.

Convenções (problema de minimização):
- o tableau é B^-1 [A | S], uma linha por variável básica e a linha Z com os custos reduzidos;
- restrições >= têm coluna excesso_<nome>, restrições <= têm coluna folga_<nome>;
- preco_sombra é dZ/db_i no sentido em que a restrição foi escrita no arquivo.
"""

from __future__ import annotations

import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from RUFAS.biophysical.animal.ration.simplex_optimizer import (  # noqa: E402
    GREATER_EQUAL,
    LESS_EQUAL,
    LinearProgram,
    LinearProgramSolver,
    SimplexResult,
)

SIMPLEX_DIR = REPO / "output" / "simplex"
COMPOSICAO = REPO / "GO_campo_dourado" / "inputs_resolvidos" / "NASEM_Comp.csv"
EPS = 1e-9

FEED_NAMES = {
    "51": "Silagem de milho",
    "110": "Silagem de leguminosa",
    "44": "Milho grão moído",
    "167": "Farelo de soja",
    "58": "Caroço de algodão",
    "301": "Núcleo mineral",
}

CONSTRAINT_NAMES = {
    "NASEM_net_energy_constraint": "Energia líquida (Mcal/dia)",
    "protein_constraint_lower": "Proteína metabolizável mínima (g/dia)",
    "protein_constraint_upper": "Proteína metabolizável máxima (g/dia)",
    "calcium_constraint": "Cálcio absorvível (g/dia)",
    "phosphorus_constraint": "Fósforo absorvível (g/dia)",
    "NDF_constraint_lower": "FDN >= 25 % da MS",
    "NDF_constraint_upper": "FDN <= 45 % da MS",
    "forage_NDF_constraint": "FDN de forragem >= 15 % da MS",
    "fat_constraint": "Gordura <= 7 % da MS",
    "DMI_constraint_lower": "Consumo de MS mínimo (kg/dia)",
    "DMI_constraint_upper": "Consumo de MS máximo (kg/dia)",
}

FILE_PATTERN = re.compile(r"modelo_(?P<pen>.+)_dia(?P<day>\d+)_tent(?P<attempt>\d+)\.csv")


def feed_label(feed_id: str) -> str:
    name = FEED_NAMES.get(str(feed_id))
    return f"{name} ({feed_id})" if name else f"Alimento {feed_id}"


def constraint_label(name: str) -> str:
    if name.startswith("max_"):
        return f"Limite superior de {feed_label(name.removeprefix('max_'))}"
    return CONSTRAINT_NAMES.get(name, name)


@dataclass(frozen=True)
class Formulacao:
    pen: str
    day: int
    attempt: int
    modelo: Path
    tableau: Path

    @property
    def label(self) -> str:
        return f"{self.pen} - dia {self.day} - tentativa {self.attempt}"


def listar_fazendas(base: Path = SIMPLEX_DIR) -> list[Path]:
    if not base.is_dir():
        return []
    return sorted(p for p in base.iterdir() if p.is_dir() and any(p.glob("modelo_*.csv")))


def listar_formulacoes(pasta: Path) -> list[Formulacao]:
    found = []
    for modelo in pasta.glob("modelo_*.csv"):
        match = FILE_PATTERN.match(modelo.name)
        if not match:
            continue
        tableau = modelo.with_name(modelo.name.replace("modelo_", "tableau_", 1))
        if tableau.exists():
            found.append(Formulacao(match["pen"], int(match["day"]), int(match["attempt"]), modelo, tableau))
    return sorted(found, key=lambda f: (f.pen, f.day, f.attempt))


@dataclass
class Modelo:
    """Conteúdo de um modelo_*.csv."""

    restricoes: pd.DataFrame  # índice: restrição; colunas: tipo, <feeds>, sinal, b, valor, folga, preco_sombra
    variaveis: pd.DataFrame  # índice: custo_c, x_otimo, custo_reduzido, limite_inferior, limite_superior
    resumo: dict[str, str]
    feeds: list[str]

    @property
    def funcionais(self) -> pd.DataFrame:
        return self.restricoes[self.restricoes["tipo"] != "limite_superior"]

    @property
    def Z(self) -> float:
        return float(self.resumo["Z"])

    def lp(
        self,
        c: np.ndarray | None = None,
        b: np.ndarray | None = None,
        lower: np.ndarray | None = None,
        upper: np.ndarray | None = None,
    ) -> LinearProgram:
        """O PL do arquivo, com custos, lados direitos ou limites das variáveis trocados."""
        rows = self.funcionais
        var = self.variaveis.loc[:, self.feeds]
        return LinearProgram(
            c=var.loc["custo_c"].to_numpy(dtype=float) if c is None else np.asarray(c, float),
            A=rows[self.feeds].to_numpy(dtype=float),
            b=rows["b"].to_numpy(dtype=float) if b is None else np.asarray(b, float),
            senses=list(rows["sinal"]),
            row_names=list(rows.index),
            col_names=list(self.feeds),
            lower=var.loc["limite_inferior"].to_numpy(dtype=float) if lower is None else np.asarray(lower, float),
            upper=var.loc["limite_superior"].to_numpy(dtype=float) if upper is None else np.asarray(upper, float),
            row_kinds=list(rows["tipo"]),
        )


def _blocos(path: Path) -> list[list[list[str]]]:
    """Separa o CSV nos blocos delimitados por linha em branco."""
    blocos, atual = [], []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle):
            if any(cell.strip() for cell in row):
                atual.append(row)
            elif atual:
                blocos.append(atual)
                atual = []
    if atual:
        blocos.append(atual)
    return blocos


def ler_modelo(path: Path) -> Modelo:
    blocos = _blocos(path)

    header, *rows = blocos[0]
    feeds = header[2 : header.index("sinal")]
    restricoes = pd.DataFrame(rows, columns=header).set_index("restricao")
    for column in feeds + ["b", "valor", "folga", "preco_sombra"]:
        restricoes[column] = pd.to_numeric(restricoes[column])

    header, *rows = blocos[1]
    variaveis = pd.DataFrame([r[2:] for r in rows], index=[r[0] for r in rows], columns=feeds).apply(pd.to_numeric)

    resumo = {r[0]: r[1] for r in blocos[2][1:]}
    return Modelo(restricoes, variaveis, resumo, feeds)


def ler_tableau(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def resolver(lp: LinearProgram) -> SimplexResult:
    return LinearProgramSolver().solve(lp)


def tableau_para_dataframe(result: SimplexResult) -> pd.DataFrame:
    if result.tableau.size == 0:
        return pd.DataFrame()
    return pd.DataFrame(result.tableau, index=result.tableau_row_names, columns=result.tableau_col_names)


# ------------------------------------------------------------------------------------------------
# Análise de sensibilidade a partir do tableau ótimo
# ------------------------------------------------------------------------------------------------


def faixa_custos(tableau: pd.DataFrame, modelo: Modelo) -> pd.DataFrame:
    """
    Intervalo de cada c_j em que a base ótima (a dieta) não muda.

    - alimento não básico (no limite inferior): c_j pode subir sem limite e cair até c_j - custo reduzido;
    - alimento básico: razões entre a linha Z e a linha do alimento no tableau, nas colunas não básicas.
    """
    z = tableau.loc["Z"]
    basicas = [r for r in tableau.index if r != "Z"]
    nao_basicas = [c for c in tableau.columns if c != "RHS" and c not in basicas]
    linhas = []
    for feed in modelo.feeds:
        c_j = float(modelo.variaveis.loc["custo_c", feed])
        x_j = float(modelo.variaveis.loc["x_otimo", feed])
        reduzido = float(z[feed])
        if feed in basicas:
            row = tableau.loc[feed, nao_basicas]
            pos = [z[k] / row[k] for k in nao_basicas if row[k] > EPS]
            neg = [z[k] / row[k] for k in nao_basicas if row[k] < -EPS]
            c_min = c_j + (max(neg) if neg else -np.inf)
            c_max = c_j + (min(pos) if pos else np.inf)
            situacao = "básica"
        else:
            c_min, c_max, situacao = c_j - reduzido, np.inf, "não básica (no limite)"
        linhas.append(
            {
                "alimento": feed_label(feed),
                "feed": feed,
                "c_j (US$/kg MS)": c_j,
                "x* (kg MS)": x_j,
                "custo reduzido": reduzido,
                "situação": situacao,
                "c_j mínimo": c_min,
                "c_j máximo": c_max,
            }
        )
    return pd.DataFrame(linhas).set_index("alimento")


def _fator_exibicao(name: str, sinal: str, b: float, coeficientes: np.ndarray) -> int:
    """-1 quando a restrição foi gravada como -a.x >= -b e deve ser mostrada como a.x <= b."""
    if sinal == GREATER_EQUAL and b <= 0 and name.endswith("_upper") and np.all(coeficientes <= EPS):
        return -1
    return 1


def faixa_rhs(tableau: pd.DataFrame, modelo: Modelo) -> pd.DataFrame:
    """
    Intervalo de cada b_i em que a base ótima não muda e o preço-sombra continua valendo.

    A coluna da folga/excesso da restrição i no tableau é +-B^-1 e_i; a base continua viável
    enquanto RHS + delta * B^-1 e_i >= 0. As colunas *_exib trazem a restrição no sentido natural
    de leitura (máximos como <=).
    """
    basicas = [r for r in tableau.index if r != "Z"]
    rhs = tableau.loc[basicas, "RHS"].to_numpy(dtype=float)
    linhas = []
    for name, row in modelo.restricoes.iterrows():
        sinal, b, y = str(row["sinal"]), float(row["b"]), float(row["preco_sombra"])
        coluna = ("folga_" if sinal == LESS_EQUAL else "excesso_") + str(name)
        if coluna not in tableau.columns:
            continue
        d = tableau.loc[basicas, coluna].to_numpy(dtype=float)
        g = d if sinal == LESS_EQUAL else -d
        sobe = (-rhs[g < -EPS] / g[g < -EPS]).min() if np.any(g < -EPS) else np.inf
        desce = (-rhs[g > EPS] / g[g > EPS]).max() if np.any(g > EPS) else -np.inf
        b_min, b_max = b + desce, b + sobe
        fator = _fator_exibicao(str(name), sinal, b, row[modelo.feeds].to_numpy(dtype=float))
        linhas.append(
            {
                "restrição": constraint_label(str(name)),
                "nome": name,
                "tipo": row["tipo"],
                "b": b,
                "folga": abs(float(row["folga"])),
                "ativa": abs(float(row["folga"])) < 1e-6,
                "preço-sombra": y,
                "b mínimo": b_min,
                "b máximo": b_max,
                "fator": fator,
                "sinal_exib": sinal if fator == 1 else LESS_EQUAL,
                "b_exib": fator * b,
                "valor_exib": fator * float(row["valor"]),
                "preço-sombra_exib": fator * y,
                "b mínimo_exib": min(fator * b_min, fator * b_max),
                "b máximo_exib": max(fator * b_min, fator * b_max),
            }
        )
    return pd.DataFrame(linhas).set_index("restrição")


def dieta(modelo: Modelo, x: np.ndarray | None = None) -> pd.DataFrame:
    valores = modelo.variaveis.loc["x_otimo", modelo.feeds].to_numpy(dtype=float) if x is None else np.asarray(x, float)
    total = float(valores.sum())
    return pd.DataFrame(
        {
            "alimento": [feed_label(f) for f in modelo.feeds],
            "kg MS/vaca/dia": valores,
            "% da MS": 100.0 * valores / total if total else np.zeros_like(valores),
        }
    ).set_index("alimento")


def materia_seca(feeds: list[str], path: Path = COMPOSICAO) -> dict[str, float]:
    """% de matéria seca de cada alimento (coluna DM da tabela de composição); 100 quando não encontrado."""
    dm = {f: 100.0 for f in feeds}
    if path.exists():
        tabela = pd.read_csv(path, usecols=["rufas_id", "DM"])
        encontrados = dict(zip(tabela["rufas_id"].astype(str), tabela["DM"].astype(float)))
        dm.update({f: encontrados[f] for f in feeds if f in encontrados and encontrados[f] > 0})
    return dm
