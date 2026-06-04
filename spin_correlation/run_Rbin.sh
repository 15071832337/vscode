#!/usr/bin/env bash
set -euo pipefail

# ---- user config ----
MACRO="Analysis_Rbin.C"   # file that contains: void final_compact(...)
#INFILE="../data_pp_NoNUA/pp_data_final.root"
#INFILE="AnalysisResults_unweight_675285.root"
#INFILE="AnalysisResults_weight.root"
#INFILE="AnalysisResults_unweight_setting5.root"
INFILE="AnalysisResults_unweight_675285_MEV3.root"


TASKPATH="lambdaspincorrderived"

INTERACTIVE="false"
SAVEPNG="true"

# Controls for ME + model switches
DO_ME="true"            # true/false
ME_USE_SE_WINDOW="false" # true/false

SIGNAL_OPT="0"  # 0 = your default (2G core + CB)

#signalOpt = 0: Double Gaussian
#signalOpt = 1: Double-sided CrystalBall(DCSB)
#signalOpt = 2: Triple Gaussian
#signalOpt = 3: Triple Gaussian
#signalOpt = 4: One Gaussian + Double-sided CrystalBall(DSCB)
#signalOpt = 6: GausExpLeft/RightTail

BKG_OPT="2"     # 0 = Cheb2

#bkgOpt = 0: 3rd order Bernstein
#bkgOpt = 1: 4th order Bernstein
#bkgOpt = 2: 2nd order Chebychev
#bkgOpt = other: Exponential of 2nd order polynomial

# ---- sparse name pairs: "SE:ME" ----
pairs=(
  "hSparseLambdaLambda:hSparseLambdaLambdaMixed"
  "hSparseLambdaAntiLambda:hSparseLambdaAntiLambdaMixed"
  "hSparseAntiLambdaLambda:hSparseAntiLambdaLambdaMixed"
  "hSparseAntiLambdaAntiLambda:hSparseAntiLambdaAntiLambdaMixed"
)

ROOTBIN="root"
ROOTFLAGS="-l -b -q"

for item in "${pairs[@]}"; do
  se="${item%%:*}"
  me="${item##*:}"

  tag="${se}"  # tag based on SE name (you can change if you want)
  outdir="Output_pp_sig${SIGNAL_OPT}_bkg${BKG_OPT}_${tag}"
  mkdir -p "${outdir}"

  echo "============================================================"
  echo "Running tag=${tag}"
  echo "  SE sparse=${se}"
  echo "  ME sparse=${me}"
  echo "  outDir=${outdir}"
  echo "  doME=${DO_ME}  meUseSEWindow=${ME_USE_SE_WINDOW}"
  echo "  signalOpt=${SIGNAL_OPT}  bkgOpt=${BKG_OPT}"
  echo "============================================================"

  ${ROOTBIN} ${ROOTFLAGS} \
    "${MACRO}+(\"${INFILE}\",\"${outdir}\",\"${TASKPATH}\",\"${se}\",\"${me}\",${INTERACTIVE},${SAVEPNG},${DO_ME},${ME_USE_SE_WINDOW},${SIGNAL_OPT},${BKG_OPT})" \
    2>&1 | tee "${outdir}/log_${tag}.txt"
done

echo "All done."
