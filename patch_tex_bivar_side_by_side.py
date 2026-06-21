#!/usr/bin/env python3
"""Patch the manuscript so the two bivariate simulation figures become one full-width figure.

Usage:
  python patch_tex_bivar_side_by_side.py 'main(3).tex' main_side_by_side.tex

The patch keeps both original labels so existing \ref commands still compile.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path


COMBINED = r'''
\begin{figure}[!t]
\centering
\includegraphics[width=\textwidth]{plots/paper_bivar_efficiency_pair_Tge1000}
\caption{Bivariate linear Hawkes simulation. Left: overall normalized RMSE, with the dotted line denoting MLE parity and horizontal lines denoting population Godambe targets. Right: population asymptotic standard-error inflation by parameter, relative to the MLE benchmark.}
\label{fig:bivar-overall-efficiency}
\label{fig:bivar-godambe-parameter}
\end{figure}
'''


def replace_one(pattern: str, repl: str, text: str, label: str) -> str:
    new, n = re.subn(pattern, lambda _m: repl, text, count=1, flags=re.S)
    if n != 1:
        raise SystemExit(f"Expected to replace exactly one {label}; replaced {n}.")
    return new


def main() -> None:
    if len(sys.argv) not in {2, 3}:
        raise SystemExit("Usage: python patch_tex_bivar_side_by_side.py INPUT.tex [OUTPUT.tex]")

    inp = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) == 3 else inp.with_name(inp.stem + "_side_by_side.tex")
    text = inp.read_text()

    if "paper_bivar_efficiency_pair_Tge1000" in text:
        out.write_text(text)
        print(f"Already patched; wrote {out}")
        return

    fig1_pat = (
        r"\n\\begin\{figure\}(?:\[[^\]]*\])?\s*"
        r"\\centering\s*"
        r"\\includegraphics\[[^\]]*\]\{plots/paper_bivar_overall_efficiency_Tge1000\}\s*"
        r"\\caption\{.*?\}\s*"
        r"\\label\{fig:bivar-overall-efficiency\}\s*"
        r"\\end\{figure\}\s*"
    )
    text = replace_one(fig1_pat, "\n", text, "overall-efficiency figure")

    fig2_pat = (
        r"\n\\begin\{figure\}(?:\[[^\]]*\])?\s*"
        r"\\centering\s*"
        r"\\includegraphics\[[^\]]*\]\{plots/paper_bivar_per_parameter_godambe_se_inflation\}\s*"
        r"\\caption\{.*?\}\s*"
        r"\\label\{fig:bivar-godambe-parameter\}\s*"
        r"\\end\{figure\}\s*"
    )
    text = replace_one(fig2_pat, "\n" + COMBINED + "\n", text, "per-parameter Godambe figure")

    text = text.replace(
        "Figure~\\ref{fig:bivar-overall-efficiency} reports",
        "The left panel of Figure~\\ref{fig:bivar-overall-efficiency} reports",
    )
    text = text.replace(
        "Figure~\\ref{fig:bivar-godambe-parameter} gives",
        "The right panel of Figure~\\ref{fig:bivar-overall-efficiency} gives",
    )

    out.write_text(text)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
