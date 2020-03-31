# commands for Lasso example analysis

mmLasso4() {
  JULIA_NUM_THREADS=2 julia -O3 -p4 OutlierVsMultimodal.jl $*
}

mmLassoMC4() {
  sleep  1; mmLasso4 $* &
  sleep 30; mmLasso4 $*
  sleep 60; mmLasso4 $* &
  sleep 90; mmLasso4 $*
}

mmLassoMC4_alt5() {
  mmLassoMC4 --altmode 0.5 $*
}

mmLassoMC4_alt6() {
  mmLassoMC4 --altmode 0.6 $*
}

mmLassoMC4_alt7() {
  mmLassoMC4 --altmode 0.7 $*
}

mmLassoMC4_alt8() {
  mmLassoMC4 --altmode 0.8 $*
}

mmLassoMC4_alt9() {
  mmLassoMC4 --altmode 0.9 $*
}

mmLassoMC4AnalysisNH() {
  mmLassoMC4_alt5 --spreadNH 1.0 &
  mmLassoMC4_alt7 --spreadNH 1.0 &
  mmLassoMC4_alt9 --spreadNH 1.0
  mmLassoMC4_alt5 --spreadNH 10.0 &
  mmLassoMC4_alt7 --spreadNH 10.0 &
  mmLassoMC4_alt9 --spreadNH 10.0
  mmLassoMC4_alt5 --spreadNH 25.0 &
  mmLassoMC4_alt7 --spreadNH 25.0 &
  mmLassoMC4_alt9 --spreadNH 25.0
}
