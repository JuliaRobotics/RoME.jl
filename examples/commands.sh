# commands for Lasso example analysis

mmLasso4() {
  JULIA_NUM_THREADS=2 julia -O3 -p4 OutlierVsMultimodal.jl $*
}

mmLassoMC5() {
  sleep  1; mmLasso4 $* &
  sleep 30; mmLasso4 $* &
  sleep 60; mmLasso4 $* &
  sleep 90; mmLasso4 $* &
  sleep 120; mmLasso4 $*
}

mmLassoMC5_alt5() {
  mmLassoMC5 --altmode 0.5 $*
}

mmLassoMC5_alt6() {
  mmLassoMC5 --altmode 0.6 $*
}

mmLassoMC5_alt7() {
  mmLassoMC5 --altmode 0.7 $*
}

mmLassoMC5_alt8() {
  mmLassoMC5 --altmode 0.8 $*
}

mmLassoMC5_alt9() {
  mmLassoMC5 --altmode 0.9 $*
}

mmLassoMC5AnalysisNH() {
  sleep 000; mmLassoMC5_alt5 --spreadNH 1.0 $* &
  sleep 150; mmLassoMC5_alt7 --spreadNH 1.0 $* &
  sleep 300; mmLassoMC5_alt9 --spreadNH 1.0 $* &
  sleep 450; mmLassoMC5_alt5 --spreadNH 5.0 $* &
  sleep 600; mmLassoMC5_alt7 --spreadNH 5.0 $* &
  sleep 750; mmLassoMC5_alt9 --spreadNH 5.0 $*
  sleep 000; mmLassoMC5_alt5 --spreadNH 10.0 $* &
  sleep 150; mmLassoMC5_alt7 --spreadNH 10.0 $* &
  sleep 300; mmLassoMC5_alt9 --spreadNH 10.0 $* &
  sleep 450; mmLassoMC5_alt5 --spreadNH 25.0 $* &
  sleep 600; mmLassoMC5_alt7 --spreadNH 25.0 $* &
  sleep 750; mmLassoMC5_alt9 --spreadNH 25.0 $*
}


# find . | grep plot_x90_solve.pdf > thelist
# echo "thelist.pdf" >> thelist
# cat thelist | xargs pdfunite
# evince thelist.pdf
