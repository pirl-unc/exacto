#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
Auto-generate Quarto markdown documentation for the exacto CLI.

Builds the same argparse parser tree that ``exactolib.cli.cli_main.run``
constructs, walks the registered subparsers, and emits one ``.qmd`` per
subcommand plus an ``index.qmd`` listing them all.

Usage:
    python scripts/docs/generate_docs.py

Output is written to ``<repo_root>/docs/cli`` regardless of the current
working directory. The script never deletes anything.

Per-subcommand pages (``<subcommand>.qmd``) are regenerated on every run,
so any manual edits to those files will be lost — edit the argparse
parsers in ``python/exactolib/cli/`` instead.

The index page (``index.qmd``) is treated as user-editable: it is only
created the first time, and subsequent runs leave it alone. Edit it
freely.
"""

from __future__ import annotations

import argparse
import inspect
import re
import shlex
from pathlib import Path
from typing import Callable, Iterable

import exactolib
from exactolib.cli import cli_main


def subparser_registrars() -> list[Callable]:
    """The ``add_cli_*_arg_parser`` functions ``cli_main.run`` calls, in its order, so the generated
    docs list exactly the commands of the live CLI."""
    names = re.findall(r"sub_parsers = (add_cli_\w+)\(sub_parsers=sub_parsers\)", inspect.getsource(cli_main.run))
    return [getattr(cli_main, name) for name in names]


# Resolve paths from the script's own location so the script works regardless
# of the current working directory. Script lives at
# ``<repo_root>/scripts/docs/generate_docs.py``, so the repo root is two
# directories up.
REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_DIR = REPO_ROOT / "docs"
CLI_DOCS_DIR = DOCS_DIR / "cli"


# Concrete usage example per subcommand. Used to populate the "Example" block
# on each page alongside the auto-generated "Template" block. File names follow
# the mutant proteoform prediction pipeline (docs/pipelines): a tumor sample with
# pass-1 and pass-2 cluster-rna-reads runs in pass1/ and pass2/. Every example is
# parsed with the live parser before any page is written (check_examples).
ANNOTATION_FLAGS = """\
    --reference-gene-annotation-file gencode.gtf.gz \\
    --reference-gene-annotation-source gencode \\
    --reference-gene-annotation-assembly hg38 \\
    --reference-gene-annotation-version v45"""

EXAMPLES: dict[str, str] = {
    "annotate-vars": f"""\
exacto annotate-vars \\
    --tsv-file tumor_somatic_dna_variants.tsv \\
{ANNOTATION_FLAGS} \\
    --output-tsv-file tumor_somatic_dna_variants_annotated.tsv""",
    "build-genome-var-graph": """\
exacto build-genome-var-graph \\
    --variants-tsv-file tumor_somatic_dna_variants.tsv \\
    --fasta-file reference_genome.fasta \\
    --output-fasta-file tumor_genome_var_graph.fasta \\
    --sequence-prefix tumor""",
    "build-transcriptome-var-graph": """\
exacto build-transcriptome-var-graph \\
    --transcript-structures-tsv-file rna_variants/tumor_exacto_assembled_transcript_model_alignments.tsv \\
    --fasta-file reference_genome.fasta \\
    --output-fasta-file tumor_transcriptome_var_graph.fasta""",
    "call-germline-dna-vars": """\
exacto call-germline-dna-vars \\
    --bam-file tumor_dna.sorted.bam \\
    --bai-file tumor_dna.sorted.bam.bai \\
    --fasta-file reference_genome.fasta \\
    --preset pb \\
    --output-tsv-file tumor_germline_dna_variants.tsv""",
    "call-somatic-dna-vars": """\
exacto call-somatic-dna-vars \\
    --bam-file tumor_dna.sorted.bam \\
    --bai-file tumor_dna.sorted.bam.bai \\
    --control-bam-files normal_dna.sorted.bam \\
    --control-bai-files normal_dna.sorted.bam.bai \\
    --fasta-file reference_genome.fasta \\
    --preset pb \\
    --output-tsv-file tumor_somatic_dna_variants.tsv""",
    "call-rna-transcript-vars": f"""\
exacto call-rna-transcript-vars \\
    --bam-file stitched/tumor_exacto_stitched_transcripts.sorted.bam \\
    --reference-genome-fasta-file reference_genome.fasta \\
{ANNOTATION_FLAGS} \\
    --output-dir rna_variants/ \\
    --output-prefix tumor""",
    "call-peptide-vars": """\
exacto call-peptide-vars \\
    --proteoforms-tsv-file translation/tumor_exacto_proteoforms.tsv \\
    --reference-fasta-file reference_proteome.fasta \\
    --output-tsv-file tumor_peptide_variants.tsv \\
    --output-fasta-file tumor_peptide_variants.fasta""",
    "cluster-rna-reads": f"""\
exacto cluster-rna-reads \\
    --bam-file tumor_rna.sorted.bam \\
    --bai-file tumor_rna.sorted.bam.bai \\
    --reference-genome-fasta-file reference_genome.fasta \\
{ANNOTATION_FLAGS} \\
    --output-dir pass1/ \\
    --output-prefix tumor_pass1 \\
    --preset pb""",
    "correct-rna-reads": f"""\
exacto correct-rna-reads \\
    --bam-file tumor_rna.sorted.bam \\
    --clusters-tsv-file pass1/tumor_pass1_exacto_rna_clusters.tsv \\
    --cluster-reference-transcripts-tsv-file pass1/tumor_pass1_exacto_rna_clusters_reference_transcripts.tsv \\
    --cluster-splice-junctions-tsv-file pass1/tumor_pass1_exacto_rna_clusters_splice_junctions.tsv \\
    --cluster-variants-tsv-file pass1/tumor_pass1_exacto_rna_clusters_variants_passed.tsv \\
{ANNOTATION_FLAGS} \\
    --output-dir pass1/ \\
    --output-prefix tumor_pass1""",
    "determine-rna-consensus": """\
exacto determine-rna-consensus \\
    --tsv-file pass2/tumor_pass2_exacto_rna_clusters.tsv \\
    --fastq-file pass1/tumor_pass1_exacto_rna_corrected_reads.fastq.gz \\
    --output-dir pass2/ \\
    --output-prefix tumor_pass2""",
    "integrate-vars": f"""\
exacto integrate-vars \\
    --dna-variants-tsv-file tumor_somatic_dna_variants.tsv \\
    --rna-variants-tsv-file rna_variants/tumor_exacto_assembled_transcript_variants.tsv \\
{ANNOTATION_FLAGS} \\
    --output-tsv-file tumor_integrated_variants.tsv""",
    "quantify-rna-abundances": """\
exacto quantify-rna-abundances \\
    --clusters-tsv-file pass2/tumor_pass2_exacto_rna_clusters.tsv \\
    --cluster-reference-matches-tsv-file pass2/tumor_pass2_exacto_rna_clusters_reference_transcripts.tsv \\
    --cluster-splice-junctions-tsv-file pass2/tumor_pass2_exacto_rna_clusters_splice_junctions.tsv \\
    --cluster-variants-tsv-file pass2/tumor_pass2_exacto_rna_clusters_variants_passed.tsv \\
    --consensus-rna-reference-matches-tsv-file rna_variants/tumor_exacto_assembled_transcript_reference_transcript_matches.tsv \\
    --output-dir abundances/ \\
    --output-prefix tumor""",
    "remove-unspliced-rnas": f"""\
exacto remove-unspliced-rnas \\
    --bam-file tumor_rna.sorted.bam \\
    --bam-bai-file tumor_rna.sorted.bam.bai \\
{ANNOTATION_FLAGS} \\
    --output-bam-file tumor_rna.spliced.sorted.bam \\
    --output-bam-bai-file tumor_rna.spliced.sorted.bam.bai""",
    "stitch-reference-transcripts": f"""\
exacto stitch-reference-transcripts \\
    --bam-file pass2/tumor_pass2_exacto_rna_consensus.sorted.bam \\
    --reference-genome-fasta-file reference_genome.fasta \\
{ANNOTATION_FLAGS} \\
    --output-dir stitched/ \\
    --output-prefix tumor""",
    "translate-seqs": """\
exacto translate-seqs \\
    --fastx-file transcripts.fasta \\
    --strategy longest_orf \\
    --output-tsv-file translations.tsv \\
    --output-fasta-file translations.fasta""",
    "translate-transcripts": """\
exacto translate-transcripts \\
    --assembled-transcript-model-alignments-tsv-file rna_variants/tumor_exacto_assembled_transcript_model_alignments.tsv \\
    --assembled-transcript-variants-tsv-file rna_variants/tumor_exacto_assembled_transcript_variants.tsv \\
    --dna-variants-tsv-file tumor_somatic_dna_variants.tsv \\
    --integrated-variants-tsv-file tumor_integrated_variants.tsv \\
    --rna-consensus-tsv-file pass2/tumor_pass2_exacto_rna_consensus.tsv \\
    --stitched-transcripts-tsv-file stitched/tumor_exacto_stitched_transcripts.tsv \\
    --strategy longest_orf \\
    --output-dir translation/ \\
    --output-prefix tumor""",
}


# I/O contract per subcommand, rendered as an "At a glance" callout under
# the Description on each page. "inputs" lists file types the command
# consumes; "outputs" lists what it produces; "next_steps" links to the
# command(s) that typically run after this one in the canonical pipeline.
CLUSTER_TABLES = ("`{prefix}_exacto_rna_clusters.tsv` (read-to-cluster assignments) and "
                  "`{prefix}_exacto_rna_clusters_{summary,reference_transcripts,splice_junctions,"
                  "template_switch,variants_passed,variants_failed}.tsv`")
RNA_TRANSCRIPT_TABLES = ("`{prefix}_exacto_assembled_transcripts.tsv` and "
                         "`{prefix}_exacto_assembled_transcript_{model_alignments,variants,exons,splice_junctions,"
                         "reference_transcript_matches,nmd_predictions,filter_status}.tsv`")

IO_CONTRACT: dict[str, dict[str, str]] = {
    "annotate-vars": {
        "inputs": "`*.tsv` (DNA or RNA variant calls), `*.gtf.gz` (gene annotation)",
        "outputs": "`*.tsv` (every column and row of the input, in the same order, with eight columns appended: `position_1_genic_region`, `position_1_annotation`, `position_1_plus_1_genic_region`, `position_1_plus_1_annotation`, `position_2_minus_1_genic_region`, `position_2_minus_1_annotation`, `position_2_genic_region`, `position_2_annotation`)",
        "next_steps": "Optional — adds gene context to a variant table for review",
    },
    "build-genome-var-graph": {
        "inputs": "`*.tsv` (DNA variants from `call-germline-dna-vars` or `call-somatic-dna-vars`), `*.fasta` (reference genome)",
        "outputs": "`*.fasta` (genome variation graph sequences)",
        "next_steps": "Endpoint — feeds downstream graph-aware analyses",
    },
    "build-transcriptome-var-graph": {
        "inputs": "`*_exacto_assembled_transcript_model_alignments.tsv` (from `call-rna-transcript-vars`), `*.fasta` (reference genome)",
        "outputs": "`*.fasta` (transcriptome variation graph sequences, one record per assembled transcript)",
        "next_steps": "Endpoint — feeds downstream graph-aware analyses",
    },
    "call-germline-dna-vars": {
        "inputs": "`*.bam`, `*.bam.bai` (long-read DNA alignment carrying `cs` tags), `*.fasta` (reference genome)",
        "outputs": "`*.tsv` (one row per germline DNA variant call)",
        "next_steps": "[`build-genome-var-graph`](build-genome-var-graph.qmd), [`annotate-vars`](annotate-vars.qmd)",
    },
    "call-somatic-dna-vars": {
        "inputs": "`*.bam`, `*.bam.bai` (case), one or more control `*.bam` / `*.bam.bai`, `*.fasta` (reference genome)",
        "outputs": "`*.tsv` (one row per somatic DNA variant call)",
        "next_steps": "[`integrate-vars`](integrate-vars.qmd), [`translate-transcripts`](translate-transcripts.qmd)",
    },
    "call-rna-transcript-vars": {
        "inputs": "`*.bam` (assembled transcripts aligned with `cs` tags — in the exacto pipeline, the realigned stitched FASTA of `stitch-reference-transcripts`), `*.fasta` (reference genome), `*.gtf.gz` (gene annotation)",
        "outputs": f"{RNA_TRANSCRIPT_TABLES} under `--output-dir`",
        "next_steps": "[`integrate-vars`](integrate-vars.qmd), [`translate-transcripts`](translate-transcripts.qmd), [`quantify-rna-abundances`](quantify-rna-abundances.qmd)",
    },
    "call-peptide-vars": {
        "inputs": "`*_exacto_proteoforms.tsv` (from `translate-transcripts`), `*.fasta` (reference proteome)",
        "outputs": "`*.tsv` (mutant peptides), `*.fasta` (mutant peptide sequences)",
        "next_steps": "Endpoint — final mutant peptide output",
    },
    "cluster-rna-reads": {
        "inputs": "`*.bam`, `*.bam.bai` (spliced long-read RNA alignment carrying `cs` tags), `*.fasta` (reference genome), `*.gtf.gz` (gene annotation)",
        "outputs": f"{CLUSTER_TABLES} under `--output-dir`",
        "next_steps": "[`correct-rna-reads`](correct-rna-reads.qmd) after pass 1 on the raw reads (`--preset pb` or `ont`); [`determine-rna-consensus`](determine-rna-consensus.qmd) after pass 2 on the realigned corrected reads (`--preset corrected --allowed-variants-tsv-file` the pass-1 `_variants_passed.tsv`)",
    },
    "correct-rna-reads": {
        "inputs": "`*.bam` (the same alignment `cluster-rna-reads` ran on; reads are recovered from their primary records in their original orientation, with their base qualities, or Q60 when the BAM has none), the four pass-1 cluster TSVs (clusters, reference transcripts, splice junctions, variants passed), `*.gtf.gz` (gene annotation)",
        "outputs": "`{prefix}_exacto_rna_corrected_reads.fastq.gz` (BGZF, one record per clustered read in cluster order, read names preserved; reads in no cluster are not written). It is written as `{prefix}_exacto_rna_corrected_reads.fastq.gz.partial` and renamed when complete, and nothing is written when a clustered read is missing from the BAM",
        "next_steps": "Realign the corrected FASTQ with minimap2 (`--cs` is required), then re-run [`cluster-rna-reads`](cluster-rna-reads.qmd) with `--preset corrected` and a distinct `--output-prefix`",
    },
    "determine-rna-consensus": {
        "inputs": "`*_exacto_rna_clusters.tsv` (the pass-2 `cluster-rna-reads` run), `*.fastq[.gz]` (the reads the clusters name — the corrected FASTQ)",
        "outputs": "`{prefix}_exacto_rna_consensus.tsv`, `{prefix}_exacto_rna_consensus.fasta` (one consensus per cluster; FASTA records named by cluster ID)",
        "next_steps": "Align the consensus FASTA with minimap2 (`--cs` is required), then [`stitch-reference-transcripts`](stitch-reference-transcripts.qmd)",
    },
    "integrate-vars": {
        "inputs": "`*.tsv` (DNA variants from `call-somatic-dna-vars`, or any callset in the variant grammar), `*_exacto_assembled_transcript_variants.tsv` (from `call-rna-transcript-vars`), `*.gtf.gz`",
        "outputs": "`*.tsv` (RNA variants linked to the DNA variants that explain them)",
        "next_steps": "[`translate-transcripts`](translate-transcripts.qmd)",
    },
    "quantify-rna-abundances": {
        "inputs": "Four pass-2 `cluster-rna-reads` TSVs (clusters, reference transcripts, splice junctions, variants passed), `*_exacto_assembled_transcript_reference_transcript_matches.tsv` (from `call-rna-transcript-vars`)",
        "outputs": "`{prefix}_exacto_rna_abundances.tsv` (per cluster), `{prefix}_exacto_rna_abundances_reference_transcripts.tsv` (per reference transcript)",
        "next_steps": "Endpoint — per-cluster and per-reference-transcript abundance estimates",
    },
    "remove-unspliced-rnas": {
        "inputs": "`*.bam`, `*.bam.bai` (coordinate-sorted RNA reads or assembled transcripts), `*.gtf.gz` (gene annotation)",
        "outputs": "Filtered `*.bam`, `*.bam.bai` (new files; an output that is one of the inputs is refused)",
        "next_steps": "[`cluster-rna-reads`](cluster-rna-reads.qmd) with `--remove-unspliced-rnas no` (it applies the same filter by default)",
    },
    "stitch-reference-transcripts": {
        "inputs": "`*.bam` (aligned RNA sequences with `cs` tags — the realigned `determine-rna-consensus` FASTA), `*.fasta` (reference genome), `*.gtf.gz` (gene annotation)",
        "outputs": "`{prefix}_exacto_stitched_transcripts.tsv`, `{prefix}_exacto_stitched_transcripts.fasta` (one record per transcript, named as in the BAM)",
        "next_steps": "Realign the stitched FASTA with minimap2, run [`call-rna-transcript-vars`](call-rna-transcript-vars.qmd) on it, then pass the stitched TSV to [`translate-transcripts`](translate-transcripts.qmd) — the stitched sequences no longer match the pre-stitch coordinates",
    },
    "translate-seqs": {
        "inputs": "`*.fasta`, `*.fastq`, or a raw sequence string",
        "outputs": "`*.tsv` (translations), `*.fasta` (peptide sequences)",
        "next_steps": "Standalone utility — not part of the main proteoform pipeline",
    },
    "translate-transcripts": {
        "inputs": "`*_exacto_rna_consensus.tsv` (from `determine-rna-consensus`) **or** an assembled transcript support TSV, optionally `*_exacto_stitched_transcripts.tsv`, the `call-rna-transcript-vars` model alignments and variants TSVs, `*.tsv` (DNA variants), `*.tsv` (integrated variants)",
        "outputs": "`{prefix}_exacto_proteoforms.tsv`, `{prefix}_exacto_proteoforms.fasta`, `{prefix}_exacto_proteoform_nucleotides.tsv`",
        "next_steps": "[`call-peptide-vars`](call-peptide-vars.qmd)",
    },
}


def build_root_parser() -> argparse.ArgumentParser:
    """Mirror ``cli_main.init_arg_parser`` + ``cli_main.run`` registration."""
    parser = argparse.ArgumentParser(description="Exacto")
    parser.add_argument(
        "--version", "-v",
        action="version",
        version=f"%(prog)s version {exactolib.__version__}",
    )
    sub_parsers = parser.add_subparsers(help="Exacto sub-commands.")
    for register in subparser_registrars():
        register(sub_parsers=sub_parsers)
    return parser


def iter_subparsers(
    parser: argparse.ArgumentParser,
) -> Iterable[tuple[str, argparse.ArgumentParser, str]]:
    """Yield (name, subparser, top-level help) for every registered subcommand."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            help_by_name = {
                choice.dest: (choice.help or "")
                for choice in action._choices_actions
            }
            for name, subparser in action.choices.items():
                yield name, subparser, help_by_name.get(name, "")


def render_type(action: argparse.Action) -> str:
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return "flag"
    t = action.type
    if t is None:
        return "str"
    return getattr(t, "__name__", str(t))


def yaml_escape(text: str) -> str:
    """Make a string safe to embed inside a double-quoted YAML scalar."""
    return text.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ").strip()


_DEFAULT_PAREN_RE = re.compile(r"\s*\(default[^)]*\)", re.IGNORECASE)
_CHOICES_PAREN_RE = re.compile(r"\s*\(choices[^)]*\)", re.IGNORECASE)


def _format_default(value: object) -> str:
    """Compact default-value rendering (lists rendered as comma-joined)."""
    if value is None or value is argparse.SUPPRESS:
        return ""
    if isinstance(value, list):
        if not value:
            return ""
        return "`" + ", ".join(str(v) for v in value) + "`"
    if value == "":
        return ""
    return f"`{value}`"


def _clean_help(help_text: str | None) -> str:
    """Strip redundant '(default: ...)' / '(choices: ...)' parentheticals."""
    if not help_text:
        return ""
    text = help_text.strip().replace("\n", " ")
    text = _DEFAULT_PAREN_RE.sub("", text)
    text = _CHOICES_PAREN_RE.sub("", text)
    # Heal artifacts left by stripping the parenthetical: " ." → "." and any
    # collapsed double-spaces.
    text = re.sub(r"\s+([.,;:])", r"\1", text)
    text = re.sub(r"\s+", " ", text)
    text = text.strip().rstrip(".").rstrip()
    return text + "." if text else ""


def render_args_table(
    actions: list[argparse.Action], include_default: bool
) -> str:
    """Render argparse actions as a Quarto markdown table.

    Choices are folded into the Type column as ``str (a|b)`` so the table stays
    compact. Required-arg tables omit Default; optional-arg tables include it.
    """
    headers = ["Flag", "Type"]
    if include_default:
        headers.append("Default")
    headers.append("Description")

    rows = ["| " + " | ".join(headers) + " |"]
    rows.append("|" + "|".join([":---"] * len(headers)) + "|")

    for action in actions:
        flag = ", ".join(f"`{s}`" for s in action.option_strings) or f"`{action.dest}`"
        type_str = render_type(action)
        if action.choices:
            choices_str = "\\|".join(str(c) for c in action.choices)
            type_cell = f"`{type_str}` ({choices_str})"
        else:
            type_cell = f"`{type_str}`"
        desc_cell = _clean_help(action.help).replace("|", "\\|")

        row = [flag, type_cell]
        if include_default:
            row.append(_format_default(action.default))
        row.append(desc_cell)
        rows.append("| " + " | ".join(row) + " |")

    # Column widths sum to 100. Required tables (no Default) give Description
    # more room.
    if include_default:
        colwidths = "[28, 14, 12, 46]"
    else:
        colwidths = "[30, 16, 54]"

    rows.append("")
    rows.append(f': {{.striped .hover tbl-colwidths="{colwidths}"}}')
    return "\n".join(rows)


def render_io_contract(name: str) -> str:
    """Render the per-subcommand I/O contract as a Quarto callout block."""
    contract = IO_CONTRACT.get(name)
    if not contract:
        return ""
    return "\n".join([
        '::: {.callout-note title="At a glance"}',
        f"**Inputs:** {contract['inputs']}",
        "",
        f"**Outputs:** {contract['outputs']}",
        "",
        f"**Typical next step:** {contract['next_steps']}",
        ":::",
        "",
    ])


def render_usage(name: str, subparser: argparse.ArgumentParser) -> str:
    """Build a multi-line ``exacto <name> ...`` snippet listing required args
    (as ``--flag <placeholder>``) and optional args (as ``[--flag PLACEHOLDER]``,
    matching argparse's default usage style)."""
    actions = [
        a for a in subparser._actions
        if not isinstance(a, argparse._HelpAction)
    ]
    if not actions:
        return f"exacto {name}"

    parts = [f"exacto {name}"]

    # Required first, then optional — same convention argparse uses.
    for action in actions:
        if not action.required:
            continue
        flag = action.option_strings[0] if action.option_strings else action.dest
        if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            placeholder = ""
        elif action.choices:
            placeholder = "<" + "|".join(str(c) for c in action.choices) + ">"
        else:
            placeholder = f"<{action.dest}>"
        parts.append(f"    {flag} {placeholder}".rstrip())

    for action in actions:
        if action.required:
            continue
        flag = action.option_strings[0] if action.option_strings else action.dest
        if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            chunk = f"[{flag}]"
        else:
            placeholder = action.dest.upper()
            if action.nargs in ("*", "+"):
                placeholder = f"{placeholder} [{placeholder} ...]"
            chunk = f"[{flag} {placeholder}]"
        parts.append(f"    {chunk}")

    return " \\\n".join(parts)


def render_subcommand(
    name: str, subparser: argparse.ArgumentParser, help_text: str
) -> str:
    description = (subparser.description or help_text or "").strip()
    lines = [
        "---",
        f'title: "{name}"',
        f'description: "{yaml_escape(description)}"',
        "---",
        "",
        "## Usage",
        "",
        "Template:",
        "",
        "```bash",
        render_usage(name, subparser),
        "```",
        "",
    ]

    example = EXAMPLES.get(name)
    if example:
        lines += [
            "Example:",
            "",
            "```bash",
            example,
            "```",
            "",
        ]

    if description:
        lines += ["## Description", "", description, ""]

    contract_block = render_io_contract(name)
    if contract_block:
        lines.append(contract_block)

    for group in subparser._action_groups:
        actions = [
            a for a in group._group_actions
            if not isinstance(a, argparse._HelpAction)
        ]
        if not actions:
            continue
        title = (group.title or "Arguments").strip()
        title = title[:1].upper() + title[1:]
        all_required = all(a.required for a in actions)
        lines += [
            f"## {title}",
            "",
            render_args_table(actions, include_default=not all_required),
            "",
        ]

    return "\n".join(lines).rstrip() + "\n"


def render_index(
    subcommands: list[tuple[str, argparse.ArgumentParser, str]],
) -> str:
    lines = [
        "---",
        'title: "Mutant Proteoform Prediction Pipeline"',
        "---",
        "",
        f"`exacto` version {exactolib.__version__}.",
        "",
        "## Subcommands",
        "",
    ]
    for name, _subparser, help_text in sorted(subcommands, key=lambda t: t[0]):
        suffix = f" — {help_text.strip()}" if help_text.strip() else ""
        lines.append(f"- [`{name}`]({name}.qmd){suffix}")
    return "\n".join(lines) + "\n"


def check_examples(parser: argparse.ArgumentParser) -> None:
    """Parse every example with the live parser, so an example that names a removed flag or misses a
    required one stops the run instead of reaching the docs."""
    for name, example in EXAMPLES.items():
        argv = shlex.split(example.replace("\\\n", " "))
        assert argv[0] == "exacto" and argv[1] == name, f"example for {name} runs {argv[:2]}"
        try:
            parser.parse_args(argv[1:])
        except SystemExit:
            raise SystemExit(f"The example for {name} does not parse with the current CLI (see the error above).")


def main() -> int:
    print("Generating Quarto CLI documentation...")
    CLI_DOCS_DIR.mkdir(parents=True, exist_ok=True)

    parser = build_root_parser()
    subcommands = list(iter_subparsers(parser))
    missing = sorted({name for name, _, _ in subcommands} - set(EXAMPLES))
    assert not missing, f"no example for {missing}"
    check_examples(parser)

    for name, subparser, help_text in subcommands:
        out_path = CLI_DOCS_DIR / f"{name}.qmd"
        out_path.write_text(render_subcommand(name, subparser, help_text))
        print(f"  wrote {out_path.relative_to(REPO_ROOT)}")

    index_path = CLI_DOCS_DIR / "index.qmd"
    if index_path.exists():
        print(f"  kept {index_path.relative_to(REPO_ROOT)} (user-editable; not overwritten)")
    else:
        index_path.write_text(render_index(subcommands))
        print(f"  wrote {index_path.relative_to(REPO_ROOT)} (initial scaffold; safe to edit)")

    print(f"Done! {len(subcommands)} subcommand pages generated.")
    print(f"  Preview with: cd {DOCS_DIR.relative_to(REPO_ROOT)} && quarto preview")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
