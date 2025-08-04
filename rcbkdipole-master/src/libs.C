void libs() {

  gROOT->ProcessLine(".L tools.cpp++");
  gROOT->ProcessLine(".L datafile.cpp++");
  gSystem->AddIncludePath("-I/opt/homebrew/opt/gsl/include");
  gSystem->Load("/opt/homebrew/opt/gsl/lib/libgsl.28.dylib");
  gROOT->ProcessLine(".L interpolation.cpp++");
  gROOT->ProcessLine(".L amplitudelib.cpp++");
  
  gROOT->ProcessLine(".L /Users/suyoupeng/Downloads/rcbkdipole-master/src/qcd.cpp++");

  gROOT->ProcessLine(".L /Users/suyoupeng/Downloads/amplitudelib-2/amplitudelib/qcd.cpp++");
  
  gROOT->ProcessLine(".L /Users/suyoupeng/Downloads/rcbkdipole-master/src/qcd.hpp++");
  gROOT->ProcessLine(".L /Users/suyoupeng/Downloads/rcbkdipole-master/src/config.hpp++");
  gROOT->ProcessLine(".L /Users/suyoupeng/Downloads/rcbkdipole-master/src/single_inclusive.hpp++");
  gROOT->ProcessLine(".L /Users/suyoupeng/Downloads/rcbkdipole-master/src/pdf.cpp++");
  gROOT->ProcessLine(".L /Users/suyoupeng/Downloads/rcbkdipole-master/src/fragmentation.hpp++");

 

  
 
  gROOT->ProcessLine(".L /Users/suyoupeng/Downloads/rcbkdipole-master/src/kkp.cpp++");
  
  //gSystem->Load("/Users/suyoupeng/Downloads/rcbkdipole-master/src/single_inclusive_hpp.so");

  //gSystem->Load("/Users/suyoupeng/Downloads/rcbkdipole-master/src/libkkp.dylib");

  //gfortran -c -O2 -fPIC -arch arm64 fragmentation_kkp.f -o kkp.o
  //gfortran -dynamiclib -arch arm64 kkp.o -o libkkp.dylib

  //  gSystem->ComileMacro();
}


//

