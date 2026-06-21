#!/bin/bash
set -euo pipefail
RESULTS_DIR="${1:-bivar_results}"
POP_DIR="${2:-bivar_population}"
FIG_DIR="${3:-bivar_figures}"
PAPER_FIG_DIR="${4:-bivar_figures_paper}"

echo "Estimator result directory: ${RESULTS_DIR}"
echo "  raw arrays:       $(find "${RESULTS_DIR}" -maxdepth 1 -name 'raw_array*.csv' 2>/dev/null | wc -l | tr -d ' ')"
echo "  manifests:        $(find "${RESULTS_DIR}" -maxdepth 1 -name 'manifest_array*.csv' 2>/dev/null | wc -l | tr -d ' ')"
echo "  failure files:    $(find "${RESULTS_DIR}" -maxdepth 1 -name 'failures_array*.csv' 2>/dev/null | wc -l | tr -d ' ')"
if compgen -G "${RESULTS_DIR}/failures_array*.csv" >/dev/null; then
  echo "  non-empty failures:"
  for f in ${RESULTS_DIR}/failures_array*.csv; do
    [ -s "$f" ] && echo "    $f ($(($(wc -l < "$f")-1)) rows)"
  done
fi

echo "Population result directory: ${POP_DIR}"
echo "  matrix arrays:    $(find "${POP_DIR}" -maxdepth 1 -name 'population_matrices_array*.csv' 2>/dev/null | wc -l | tr -d ' ')"
echo "  summary arrays:   $(find "${POP_DIR}" -maxdepth 1 -name 'population_summary_array*.csv' 2>/dev/null | wc -l | tr -d ' ')"
echo "  failure files:    $(find "${POP_DIR}" -maxdepth 1 -name 'population_failures_array*.csv' 2>/dev/null | wc -l | tr -d ' ')"
if compgen -G "${POP_DIR}/population_failures_array*.csv" >/dev/null; then
  echo "  non-empty population failures:"
  for f in ${POP_DIR}/population_failures_array*.csv; do
    [ -s "$f" ] && echo "    $f ($(($(wc -l < "$f")-1)) rows)"
  done
fi

echo "Combined outputs:"
for f in \
  "${RESULTS_DIR}/bivar_design.csv" \
  "${RESULTS_DIR}/bivar_raw_combined.csv" \
  "${RESULTS_DIR}/population_godambe_summary.csv" \
  "${RESULTS_DIR}/population_parameter_inflation.csv" \
  "${RESULTS_DIR}/population_parameter_inflation_long.csv" \
  "${RESULTS_DIR}/population_relative_eigenvalues.csv" \
  "${RESULTS_DIR}/bivar_standardized_summary.csv" \
  "${RESULTS_DIR}/bivar_asymptotic_ci_by_parameter.csv" \
  "${RESULTS_DIR}/bivar_asymptotic_ci_overall.csv" \
  "${RESULTS_DIR}/bivar_lan_quadratic_quantiles.csv" \
  "${RESULTS_DIR}/bivar_rate_slopes.csv"; do
  [ -e "$f" ] && echo "  present: $f" || echo "  missing: $f"
done

echo "Diagnostic PDF figure outputs:"
if [ -d "${FIG_DIR}" ]; then
  find "${FIG_DIR}" -maxdepth 1 -name '*.pdf' -print | sort
else
  echo "  ${FIG_DIR} not found"
fi

echo "Paper PDF figure outputs:"
if [ -d "${PAPER_FIG_DIR}" ]; then
  find "${PAPER_FIG_DIR}" -maxdepth 1 -name '*.pdf' -print | sort
else
  echo "  ${PAPER_FIG_DIR} not found"
fi

echo "Paper PNG figure outputs:"
if [ -d "${PAPER_FIG_DIR}" ]; then
  find "${PAPER_FIG_DIR}" -maxdepth 1 -name '*.png' -print | sort
else
  echo "  ${PAPER_FIG_DIR} not found"
fi

pair="${PAPER_FIG_DIR}/paper_bivar_efficiency_pair_Tge1000.pdf"
pilot_pair="${PAPER_FIG_DIR}/paper_bivar_pilot_efficiency_pair_Tge300.pdf"
if [ -e "$pair" ]; then
  echo "Manuscript side-by-side figure: present: $pair"
elif [ -e "$pilot_pair" ]; then
  echo "Pilot side-by-side figure: present: $pilot_pair"
else
  echo "Side-by-side figure: missing expected full or pilot pair in ${PAPER_FIG_DIR}"
fi
