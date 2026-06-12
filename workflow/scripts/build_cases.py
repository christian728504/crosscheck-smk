"""Stage E: build the sample-swap "cases" table and per-cluster Bokeh graphs.

From the aggregated crosscheck results, keep only WGS-on-left cross-modality
comparisons, then retain every UNEXPECTED_* pair plus the EXPECTED_MATCH pairs that
share a sample with an unexpected one (within the same comparison) -- so each anomaly
is shown next to the correct relationship it conflicts with. The kept pairs form a
graph (samples = nodes, comparisons = edges); each connected component is one "case".

Outputs:
  - results/cases.tsv          : the kept pairs, annotated with SUBGRAPH_ID.
  - results/cases/case_{n}.html: one self-contained Bokeh plot per component,
                                 numbered largest-first.

Node labels are participant_ids; the MOHD accession (entity:dataset_id) is in the
node hover tooltip. Modality and participant come from the manifest, not sample-name
parsing. Run via Snakemake's ``script:`` directive (``snakemake`` object injected).
"""

from pathlib import Path

import networkx as nx
import polars as pl
from bokeh.io import save
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from bokeh.resources import INLINE

# Edge color by crosscheck result; node color by modality (manifest labels).
RESULT_COLOR = {
    "EXPECTED_MATCH": "#2ca02c",
    "UNEXPECTED_MISMATCH": "#d62728",
    "UNEXPECTED_MATCH": "#ff7f0e",
}
MOD_COLOR = {
    "WGS": "#1f77b4",
    "WGBS": "#9467bd",
    "RNA-seq": "#17becf",
    "ATAC-seq": "#bcbd22",
}


def layout(graph):
    """graphviz `sfdp` positions (even, symmetric); spring_layout fallback."""
    try:
        return nx.nx_pydot.graphviz_layout(graph, prog="sfdp")
    except Exception as exc:  # noqa: BLE001 - any graphviz/pydot failure
        print(f"graphviz sfdp unavailable ({exc}); falling back to spring_layout")
        return nx.spring_layout(graph, seed=0)


def plot_case(graph, component, pos, participant_of, mod_of, path):
    """Render one connected component to a standalone Bokeh HTML file."""
    sub = graph.subgraph(component)

    plot = figure(
        height=850,
        sizing_mode="stretch_width",
        hidpi=True,
        x_axis_location=None,
        y_axis_location=None,
        title=f"{path.stem}  ({sub.number_of_nodes()} samples, {sub.number_of_edges()} comparisons)",
        tools="pan,wheel_zoom,box_zoom,reset,save",
        background_fill_color="#fafafa",
    )
    plot.grid.grid_line_color = None

    # Edges, grouped by result so each gets its own color + legend entry.
    edge_renderers = []
    for result, color in RESULT_COLOR.items():
        xs, ys, lefts, rights, lods, comps = [], [], [], [], [], []
        for u, v, d in sub.edges(data=True):
            if d["result"] != result:
                continue
            xs.append([pos[u][0], pos[v][0]])
            ys.append([pos[u][1], pos[v][1]])
            lefts.append(u)
            rights.append(v)
            lods.append(d["lod"])
            comps.append(d["comparison"])
        if not xs:
            continue
        src = ColumnDataSource(
            dict(xs=xs, ys=ys, left=lefts, right=rights, lod=lods,
                 comparison=comps, result=[result] * len(xs))
        )
        edge_renderers.append(
            plot.multi_line("xs", "ys", source=src, line_color=color,
                            line_alpha=0.85, line_width=2, legend_label=result)
        )

    # Nodes, grouped by modality.
    node_renderers = []
    for modality, color in MOD_COLOR.items():
        nodes = [n for n in sub.nodes if mod_of.get(n) == modality]
        if not nodes:
            continue
        src = ColumnDataSource(
            dict(
                x=[pos[n][0] for n in nodes],
                y=[pos[n][1] for n in nodes],
                accession=nodes,
                participant=[participant_of.get(n) for n in nodes],
                modality=[modality] * len(nodes),
                degree=[sub.degree[n] for n in nodes],
            )
        )
        node_renderers.append(
            plot.scatter("x", "y", source=src, size=18, fill_color=color,
                         line_color="#222", legend_label=modality)
        )
        # Label each node with its participant_id (offset above the marker in px).
        plot.text("x", "y", text="participant", source=src, y_offset=-12,
                  text_align="center", text_baseline="bottom", text_font_size="10pt")

    if node_renderers:
        plot.add_tools(HoverTool(renderers=node_renderers, tooltips=[
            ("accession", "@accession"),
            ("participant", "@participant"),
            ("modality", "@modality"),
            ("edges", "@degree"),
        ]))
    if edge_renderers:
        plot.add_tools(HoverTool(renderers=edge_renderers, line_policy="interp", tooltips=[
            ("pair", "@left  ↔  @right"),
            ("result", "@result"),
            ("LOD", "@lod{0.00}"),
            ("comparison", "@comparison"),
        ]))

    plot.legend.location = "top_right"
    plot.legend.click_policy = "hide"
    save(plot, filename=str(path), resources=INLINE, title=path.stem)


def main():
    combined = pl.read_parquet(snakemake.input.combined)  # noqa: F821
    manifest = pl.read_parquet(snakemake.input.manifest)  # noqa: F821

    mod_of = dict(zip(
        manifest.get_column("entity:dataset_id").to_list(),
        manifest.get_column("data_modality").to_list(),
    ))
    participant_of = dict(zip(
        manifest.get_column("entity:dataset_id").to_list(),
        manifest.get_column("participant_id").to_list(),
    ))

    # WGS-on-left, cross-modality comparisons only (drops WGS x WGS, same-modality,
    # and the duplicate WGS-on-right orientation). Modality from the manifest join.
    crosscheck_df = (
        combined.select("LEFT_GROUP_VALUE", "RIGHT_GROUP_VALUE", "RESULT",
                        "LOD_SCORE", "COMPARISON")
        .with_columns(
            pl.col("LEFT_GROUP_VALUE").replace_strict(mod_of, default=None).alias("LEFT_MOD"),
            pl.col("RIGHT_GROUP_VALUE").replace_strict(mod_of, default=None).alias("RIGHT_MOD"),
        )
        .filter((pl.col("LEFT_MOD") == "WGS") & (pl.col("LEFT_MOD") != pl.col("RIGHT_MOD")))
    )

    # Keep UNEXPECTED_* plus the EXPECTED_MATCH that share a node with an unexpected
    # edge within the same comparison (left matched against unexpected lefts, right
    # against unexpected rights). Done with semi-joins on (COMPARISON, node) -- the
    # window+implode form is rejected by recent polars.
    unexpected = crosscheck_df.filter(pl.col("RESULT").str.starts_with("UNEXPECTED_"))
    cols = crosscheck_df.columns
    unexp_left = (
        unexpected.select("COMPARISON", "LEFT_GROUP_VALUE").unique()
        .with_columns(pl.lit(True).alias("_lflag"))
    )
    unexp_right = (
        unexpected.select("COMPARISON", "RIGHT_GROUP_VALUE").unique()
        .with_columns(pl.lit(True).alias("_rflag"))
    )
    expected_keep = (
        crosscheck_df.filter(pl.col("RESULT") == "EXPECTED_MATCH")
        .join(unexp_left, on=["COMPARISON", "LEFT_GROUP_VALUE"], how="left")
        .join(unexp_right, on=["COMPARISON", "RIGHT_GROUP_VALUE"], how="left")
        .filter(pl.col("_lflag").is_not_null() | pl.col("_rflag").is_not_null())
        .select(cols)
    )
    stacked_df = pl.concat([unexpected.select(cols), expected_keep], how="vertical")

    graph = nx.Graph()
    for r in stacked_df.iter_rows(named=True):
        graph.add_edge(r["LEFT_GROUP_VALUE"], r["RIGHT_GROUP_VALUE"],
                       result=r["RESULT"], lod=r["LOD_SCORE"], comparison=r["COMPARISON"])

    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    subgraph_id = {n: i for i, comp in enumerate(components) for n in comp}

    # Annotate the table with participant ids and the SUBGRAPH_ID (by left node).
    stacked_df = stacked_df.with_columns(
        pl.col("LEFT_GROUP_VALUE").replace_strict(participant_of, default=None).alias("LEFT_PARTICIPANT"),
        pl.col("RIGHT_GROUP_VALUE").replace_strict(participant_of, default=None).alias("RIGHT_PARTICIPANT"),
        pl.col("LEFT_GROUP_VALUE")
        .replace_strict(subgraph_id, default=None, return_dtype=pl.Int32)
        .alias("SUBGRAPH_ID"),
    ).select(
        "SUBGRAPH_ID", "COMPARISON", "RESULT", "LOD_SCORE",
        "LEFT_GROUP_VALUE", "LEFT_PARTICIPANT", "LEFT_MOD",
        "RIGHT_GROUP_VALUE", "RIGHT_PARTICIPANT", "RIGHT_MOD",
    ).sort("SUBGRAPH_ID", "RESULT")

    out_dir = Path(snakemake.output.plots)  # noqa: F821
    out_dir.mkdir(parents=True, exist_ok=True)
    stacked_df.write_csv(snakemake.output.tsv, separator="\t")  # noqa: F821

    print(f"{stacked_df.height} kept pairs across {len(components)} cases")
    for i, component in enumerate(components):
        pos = layout(graph.subgraph(component))
        plot_case(graph, component, pos, participant_of, mod_of, out_dir / f"case_{i}.html")


if __name__ == "__main__":
    main()
