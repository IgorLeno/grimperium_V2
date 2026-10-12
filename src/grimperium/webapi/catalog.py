"""What a routine can be built from: properties, input modes, programs, outputs.

Every default and every choice list is read from the Semi-Imperium code
(settings dataclasses and the CREST/MOPAC runners), so the dashboard can
never offer a value the backend would reject. Executables are probed on
PATH at request time; nothing is launched.
"""

from __future__ import annotations

from grimperium.webapi.models import (
    RoutineCatalog,
    RoutineParameter,
    RoutineProgram,
    RoutineProperty,
)
from semi_imperium.conformers.crest_runner import (
    _METHOD_FLAGS,
    _PREOPTIMIZERS,
    _QUICK_FLAGS,
)
from semi_imperium.conformers.initial_structure import RDKitInitialStructure
from semi_imperium.domain.configuration import (
    ConformerSearchSettings,
    ConformerSelectionSettings,
    FoldingFilterSettings,
    SemiempiricalSettings,
    VerificationSettings,
)
from semi_imperium.domain.enums import ConformerSelectionStrategy, VerificationPolicy
from semi_imperium.mopac.models import SUPPORTED_HAMILTONIANS
from semi_imperium.settings import PROPERTY_ID, RuntimeSettings

P = RoutineParameter


def routine_catalog() -> RoutineCatalog:
    search = ConformerSearchSettings()
    selection = ConformerSelectionSettings()
    folding = FoldingFilterSettings()
    mopac = SemiempiricalSettings()
    verification = VerificationSettings()
    runtime = RuntimeSettings()
    rdkit = RDKitInitialStructure()
    tools = {tool.name: tool for tool in runtime.readiness()}

    properties = [
        RoutineProperty(
            key=PROPERTY_ID,
            label="Entalpia de formação (ΔHf° 298 K)",
            unit="kcal/mol",
            description="Energia para formar a molécula a partir dos elementos; "
            "única propriedade calculada hoje.",
            available=True,
        )
    ]

    input_modes = [
        RoutineProgram(
            key="smiles",
            label="SMILES",
            role="input",
            description="Uma molécula por linha em notação SMILES (ex.: CCO).",
            executable=None,
            available=None,
            parameters=[],
        ),
        RoutineProgram(
            key="name",
            label="Nome da molécula (PubChem)",
            role="input",
            description="O nome é resolvido para SMILES no PubChem; nomes "
            "ambíguos pedem escolha do candidato.",
            executable=None,
            available=None,
            parameters=[
                P(
                    key="on_ambiguous",
                    label="Se houver mais de um candidato",
                    kind="choice",
                    default="ask",
                    choices=["ask", "skip"],
                    help="ask: pausa para escolher; skip: ignora a molécula.",
                )
            ],
        ),
        RoutineProgram(
            key="manual",
            label="Manual (outra notação)",
            role="input",
            description="Estrutura informada à mão quando o nome não resolve; "
            "hoje o backend aceita SMILES manual.",
            executable=None,
            available=None,
            parameters=[
                P(
                    key="notation",
                    label="Notação",
                    kind="choice",
                    default="smiles",
                    choices=["smiles"],
                    help="Outras notações (InChI, XYZ) ainda não têm conversor.",
                )
            ],
        ),
    ]

    crest_tool = tools.get("CREST")
    mopac_tool = tools.get("MOPAC")
    programs = [
        RoutineProgram(
            key="rdkit",
            label="RDKit",
            role="structure",
            description="Valida o SMILES e gera uma estrutura 3D inicial "
            "(usada sozinha quando não há busca conformacional).",
            executable=None,
            available=None,
            parameters=[
                P(
                    key="random_seed",
                    label="Semente aleatória",
                    kind="int",
                    default=rdkit.random_seed,
                    help="Mesma semente, mesma geometria inicial.",
                ),
                P(
                    key="optimize",
                    label="Pré-otimizar com UFF",
                    kind="bool",
                    default=rdkit.optimize,
                ),
                P(
                    key="max_optimization_steps",
                    label="Máximo de passos UFF",
                    kind="int",
                    default=rdkit.max_optimization_steps,
                    minimum=1,
                ),
            ],
        ),
        RoutineProgram(
            key="crest",
            label="CREST",
            role="conformers",
            description="Busca as formas (conformações) de menor energia da "
            "molécula.",
            executable=runtime.crest_executable,
            available=crest_tool.available if crest_tool else None,
            parameters=[
                P(
                    key="method",
                    label="Método",
                    kind="choice",
                    default=search.method,
                    choices=list(_METHOD_FLAGS),
                ),
                P(
                    key="quick_mode",
                    label="Modo rápido",
                    kind="choice",
                    default=search.quick_mode,
                    choices=list(_QUICK_FLAGS),
                ),
                P(
                    key="opt_level",
                    label="Nível de otimização",
                    kind="int",
                    default=search.opt_level,
                ),
                P(
                    key="rmsd_threshold",
                    label="Limiar RMSD",
                    kind="float",
                    default=search.rmsd_threshold,
                    unit="Å",
                    minimum=0,
                ),
                P(
                    key="energy_window_kcal_mol",
                    label="Janela de energia",
                    kind="float",
                    default=search.energy_window_kcal_mol,
                    unit="kcal/mol",
                    minimum=0,
                ),
                P(
                    key="max_conformers",
                    label="Máximo de conformações",
                    kind="int",
                    default=search.max_conformers,
                    minimum=1,
                ),
                P(key="nci", label="Modo NCI", kind="bool", default=search.nci),
                P(key="use_v3", label="CREST v3", kind="bool", default=search.use_v3),
                P(
                    key="preoptimizer",
                    label="Pré-otimizador",
                    kind="choice",
                    default=search.preoptimizer,
                    choices=list(_PREOPTIMIZERS),
                ),
                P(
                    key="crest_threads",
                    label="Threads",
                    kind="int",
                    default=runtime.crest_threads,
                    minimum=1,
                    help="Só afeta o tempo, não o resultado.",
                ),
                P(
                    key="crest_timeout_seconds",
                    label="Tempo limite",
                    kind="float",
                    default=runtime.crest_timeout_seconds,
                    unit="s",
                    minimum=1,
                ),
            ],
        ),
        RoutineProgram(
            key="confpass",
            label="Seleção de conformações (CONFPASS)",
            role="selection",
            description="Escolhe quais conformações seguem para o MOPAC.",
            executable=None,
            available=None,
            experimental=True,
            parameters=[
                P(
                    key="strategy",
                    label="Estratégia",
                    kind="choice",
                    default=selection.strategy,
                    choices=[item.value for item in ConformerSelectionStrategy],
                    help="confpass_prioritization é experimental.",
                ),
                P(
                    key="top_n",
                    label="Quantas conformações manter",
                    kind="int",
                    default=selection.top_n,
                    minimum=1,
                ),
                P(
                    key="energy_window_kcal_mol",
                    label="Janela extra de energia",
                    kind="float",
                    default=selection.energy_window_kcal_mol,
                    unit="kcal/mol",
                    help="Vazio: sem limite extra.",
                ),
                P(
                    key="folding_filter",
                    label="Filtro de dobramento",
                    kind="bool",
                    default=folding.enabled,
                    help="Descarta conformações excessivamente dobradas.",
                ),
            ],
        ),
        RoutineProgram(
            key="mopac",
            label="MOPAC",
            role="energy",
            description="Otimiza a geometria e calcula a entalpia semiempírica.",
            executable=runtime.mopac_executable,
            available=mopac_tool.available if mopac_tool else None,
            parameters=[
                P(
                    key="hamiltonian",
                    label="Hamiltoniano",
                    kind="choice",
                    default=mopac.hamiltonian,
                    choices=list(SUPPORTED_HAMILTONIANS),
                ),
                P(
                    key="keywords",
                    label="Palavras-chave",
                    kind="text",
                    default=" ".join(mopac.keywords),
                ),
                P(
                    key="scf_convergence",
                    label="Convergência SCF",
                    kind="float",
                    default=mopac.scf_convergence,
                ),
                P(
                    key="solvent",
                    label="Meio",
                    kind="choice",
                    default=mopac.solvent,
                    choices=[mopac.solvent],
                ),
                P(
                    key="verification",
                    label="Verificar mínimo (frequências)",
                    kind="choice",
                    default=verification.policy.value,
                    choices=[item.value for item in VerificationPolicy],
                ),
                P(
                    key="mopac_timeout_seconds",
                    label="Tempo limite",
                    kind="float",
                    default=runtime.mopac_timeout_seconds,
                    unit="s",
                    minimum=1,
                ),
            ],
        ),
    ]

    outputs = [
        RoutineProgram(
            key="spreadsheet",
            label="Planilha",
            role="output",
            description="Uma linha por molécula com valores e proveniência.",
            executable=None,
            available=None,
            parameters=[
                P(
                    key="format",
                    label="Formato",
                    kind="choice",
                    default="csv",
                    choices=["csv", "xlsx"],
                ),
                P(
                    key="include_descriptors",
                    label="Incluir descritores MOPAC/RDKit",
                    kind="bool",
                    default=True,
                ),
            ],
        ),
        RoutineProgram(
            key="value",
            label="Valor na tela",
            role="output",
            description="Mostra o ΔHf de cada molécula ao final, sem gerar arquivo.",
            executable=None,
            available=None,
            parameters=[],
        ),
    ]

    return RoutineCatalog(
        dispatch_enabled=False,
        properties=properties,
        input_modes=input_modes,
        programs=programs,
        outputs=outputs,
    )
