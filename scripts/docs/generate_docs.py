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
    python scripts/docs/generate_cli_docs.py

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
import re
from pathlib import Path
from typing import Iterable

import exactolib
from exactolib.cli.cli_annotate_vars import add_cli_annotate_vars_arg_parser
from exactolib.cli.cli_build_genome_var_graph import (
    add_cli_build_genome_var_graph_arg_parser,
)
from exactolib.cli.cli_build_transcriptome_var_graph import (
    add_cli_build_transcriptome_var_graph_arg_parser,
)
from exactolib.cli.cli_call_germline_dna_vars import (
    add_cli_call_germline_dna_vars_arg_parser,
)
from exactolib.cli.cli_call_somatic_dna_vars import (
    add_cli_call_somatic_dna_vars_arg_parser,
)
from exactolib.cli.cli_call_peptide_vars import add_cli_call_peptide_vars_arg_parser
from exactolib.cli.cli_call_rna_vars import add_cli_call_rna_vars_arg_parser
from exactolib.cli.cli_integrate_vars import add_cli_integrate_vars_arg_parser
from exactolib.cli.cli_remove_unspliced_rnas import (
    add_cli_remove_unspliced_rnas_arg_parser,
)
from exactolib.cli.cli_translate_seqs import add_cli_translate_seqs_arg_parser
from exactolib.cli.cli_translate_structs import add_cli_translate_structs_arg_parser


# Match the order in cli_main.run so generated docs reflect the live CLI.
SUBPARSER_REGISTRARS = [
    add_cli_annotate_vars_arg_parser,
    add_cli_build_genome_var_graph_arg_parser,
    add_cli_build_transcriptome_var_graph_arg_parser,
    add_cli_call_germline_dna_vars_arg_parser,
    add_cli_call_somatic_dna_vars_arg_parser,
    add_cli_call_rna_vars_arg_parser,
    add_cli_call_peptide_vars_arg_parser,
    add_cli_integrate_vars_arg_parser,
    add_cli_remove_unspliced_rnas_arg_parser,
    add_cli_translate_seqs_arg_parser,
    add_cli_translate_structs_arg_parser,
]


# Resolve paths from the script's own location so the script works regardless
# of the current working directory. Script lives at
# ``<repo_root>/scripts/docs/generate_cli_docs.py``, so the repo root is two
# directories up.
REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_DIR = REPO_ROOT / "docs"
CLI_DOCS_DIR = DOCS_DIR / "cli"


# Concrete usage example per subcommand. Used to populate the "Example" block
# on each page alongside the auto-generated "Template" block. Filenames are
# illustrative — they reflect the typical mutant-proteoform-prediction pipeline
# from the README.
EXAMPLES: dict[str, str] = {
    "annotate-vars": """\
exacto annotate-vars \\
    --tsv-file tumor_specific_dna_variants.tsv \\
    --reference-gene-annotation-file gencode.gtf.gz \\
    --reference-gene-annotation-source gencode \\
    --reference-gene-annotation-assembly hg38 \\
    --reference-gene-annotation-version v44 \\
    --output-tsv-file tumor_specific_dna_variants.annotated.tsv""",
    "build-genome-var-graph": """\
exacto build-genome-var-graph \\
    --variants-tsv-file tumor_dna_rna_variants_integrated.tsv \\
    --fasta-file reference_genome.fasta \\
    --output-dir genome_var_graph_outputs/""",
    "build-transcriptome-var-graph": """\
exacto build-transcriptome-var-graph \\
    --transcript-structures-tsv-file tumor_primary_structures.tsv \\
    --fasta-file reference_transcriptome.fasta \\
    --output-fasta-file transcriptome_var_graph.fasta""",
    "call-germline-dna-vars": """\
exacto call-germline-dna-vars \\
    --bam-file sample_dna.sorted.bam \\
    --bam-bai-file sample_dna.sorted.bam.bai \\
    --fasta-file reference_genome.fasta \\
    --preset pb \\
    --output-tsv-file sample_germline_dna_variants.tsv""",
    "call-somatic-dna-vars": """\
exacto call-somatic-dna-vars \\
    --bam-file tumor_dna.sorted.bam \\
    --bam-bai-file tumor_dna.sorted.bam.bai \\
    --control-bam-files normal_dna.sorted.bam \\
    --control-bam-bai-files normal_dna.sorted.bam.bai \\
    --fasta-file reference_genome.fasta \\
    --preset pb \\
    --output-tsv-file tumor_somatic_dna_variants.tsv""",
    "call-peptide-vars": """\
exacto call-peptide-vars \\
    --primary-structures-tsv-file tumor_primary_structures.tsv \\
    --reference-fasta-file reference_proteome.fasta \\
    --output-tsv-file tumor_peptide_variants.tsv \\
    --output-fasta-file tumor_peptide_variants.fasta""",
    "call-rna-vars": """\
exacto call-rna-vars \\
    --bam-file tumor_transcriptome_assembly.sorted.filtered.bam \\
    --bam-bai-file tumor_transcriptome_assembly.sorted.filtered.bam.bai \\
    --reference-genome-fasta-file reference_genome.fasta \\
    --reference-gene-annotation-file gencode.gtf.gz \\
    --reference-gene-annotation-source gencode \\
    --reference-gene-annotation-assembly hg38 \\
    --reference-gene-annotation-version v44 \\
    --output-dir rna_variants_outputs/ \\
    --output-prefix tumor""",
    "integrate-vars": """\
exacto integrate-vars \\
    --annotated-dna-vars-tsv-file tumor_specific_dna_variants.annotated.tsv \\
    --rna-vars-tsv-file rna_variants_outputs/tumor_exacto_rna_variant_calls.tsv \\
    --reference-gene-annotation-file gencode.gtf.gz \\
    --reference-gene-annotation-source gencode \\
    --reference-gene-annotation-assembly hg38 \\
    --reference-gene-annotation-version v44 \\
    --output-tsv-file tumor_dna_rna_variants_integrated.tsv""",
    "remove-unspliced-rnas": """\
exacto remove-unspliced-rnas \\
    --bam-file tumor_transcriptome_assembly.sorted.bam \\
    --bam-bai-file tumor_transcriptome_assembly.sorted.bam.bai \\
    --fasta-file reference_genome.fasta \\
    --reference-gene-annotation-file gencode.gtf.gz \\
    --reference-gene-annotation-source gencode \\
    --reference-gene-annotation-assembly hg38 \\
    --reference-gene-annotation-version v44 \\
    --output-bam-file tumor_transcriptome_assembly.sorted.filtered.bam \\
    --output-bam-bai-file tumor_transcriptome_assembly.sorted.filtered.bam.bai \\
    --output-fasta-file tumor_transcriptome_assembly.sorted.filtered.fasta""",
    "translate-seqs": """\
exacto translate-seqs \\
    --fasta-file tumor_transcriptome_assembly.fasta \\
    --strategy longest_orf \\
    --output-tsv-file tumor_translations.tsv \\
    --output-fasta-file tumor_translations.fasta""",
    "translate-structs": """\
exacto translate-structs \\
    --transcript-structures-tsv-file rna_variants_outputs/tumor_exacto_transcript_structures.tsv \\
    --rna-variant-calls-tsv-file rna_variants_outputs/tumor_exacto_rna_variant_calls.tsv \\
    --integrated-variants-tsv-file tumor_dna_rna_variants_integrated.tsv \\
    --strategy longest_orf \\
    --output-tsv-file tumor_primary_structures.tsv \\
    --output-fasta-file tumor_primary_structures.fasta""",
}


# I/O contract per subcommand, rendered as an "At a glance" callout under
# the Description on each page. "inputs" lists file types the command
# consumes; "outputs" lists what it produces; "next_steps" links to the
# command(s) that typically run after this one in the canonical pipeline.
IO_CONTRACT: dict[str, dict[str, str]] = {
    "annotate-vars": {
        "inputs": "`*.tsv` (DNA or RNA variant calls), `*.gtf.gz` (gene annotation)",
        "outputs": "`*.tsv` (variants with gene-level annotation)",
        "next_steps": "[`integrate-vars`](integrate-vars.qmd), [`build-genome-var-graph`](build-genome-var-graph.qmd)",
    },
    "build-genome-var-graph": {
        "inputs": "`*.tsv` (variants), `*.fasta` (reference genome)",
        "outputs": "Graph files written to `--output-dir`",
        "next_steps": "Endpoint — feeds downstream graph-aware analyses",
    },
    "build-transcriptome-var-graph": {
        "inputs": "`*.tsv` (transcript structures), `*.fasta` (reference transcriptome)",
        "outputs": "`*.fasta` (transcriptome variation graph)",
        "next_steps": "Endpoint — feeds downstream graph-aware analyses",
    },
    "call-germline-dna-vars": {
        "inputs": "`*.bam`, `*.bam.bai`, `*.fasta`",
        "outputs": "`*.tsv` (one row per germline DNA variant call)",
        "next_steps": "[`annotate-vars`](annotate-vars.qmd)",
    },
    "call-somatic-dna-vars": {
        "inputs": "`*.bam`, `*.bam.bai`, `*.fasta`, plus one or more control `*.bam` / `*.bam.bai`",
        "outputs": "`*.tsv` (one row per somatic DNA variant call)",
        "next_steps": "[`annotate-vars`](annotate-vars.qmd)",
    },
    "call-peptide-vars": {
        "inputs": "`*.tsv` (primary structures), `*.fasta` (reference proteome)",
        "outputs": "`*.tsv` (peptide variants), `*.fasta` (peptide sequences)",
        "next_steps": "Endpoint — final mutant peptide output",
    },
    "call-rna-vars": {
        "inputs": "`*.bam`, `*.bam.bai`, `*.fasta`, `*.gtf.gz`",
        "outputs": "`*.tsv` (RNA variants), `*.tsv` (transcript structures) under `--output-dir`",
        "next_steps": "[`integrate-vars`](integrate-vars.qmd), [`translate-structs`](translate-structs.qmd)",
    },
    "integrate-vars": {
        "inputs": "`*.tsv` (annotated DNA variants), `*.tsv` (RNA variants), `*.gtf.gz`",
        "outputs": "`*.tsv` (integrated DNA + RNA variants)",
        "next_steps": "[`translate-structs`](translate-structs.qmd)",
    },
    "remove-unspliced-rnas": {
        "inputs": "`*.bam`, `*.bam.bai`, `*.fasta`, `*.gtf.gz`",
        "outputs": "Filtered `*.bam`, `*.bam.bai`, `*.fasta`",
        "next_steps": "[`call-rna-vars`](call-rna-vars.qmd)",
    },
    "translate-seqs": {
        "inputs": "`*.fasta`, `*.fastq`, or a raw sequence string",
        "outputs": "`*.tsv` (translations), `*.fasta` (peptide sequences)",
        "next_steps": "Standalone utility — not part of the main proteoform pipeline",
    },
    "translate-structs": {
        "inputs": "`*.tsv` (transcript structures), `*.tsv` (RNA variants), `*.tsv` (integrated variants)",
        "outputs": "`*.tsv` (primary structures), `*.fasta` (primary structures FASTA)",
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
    for register in SUBPARSER_REGISTRARS:
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


def main() -> int:
    print("Generating Quarto CLI documentation...")
    CLI_DOCS_DIR.mkdir(parents=True, exist_ok=True)

    parser = build_root_parser()
    subcommands = list(iter_subparsers(parser))

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
