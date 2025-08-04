//Estimate dihadron correlation cross section numeric multiDim integral
//in ROOT 
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include <ctime>
#include <sstream>
#include <unistd.h>

#include "/Users/suyoupeng/Downloads/rcbkdipole-master/src/amplitudelib.hpp"

#include "TF1.h"
#include "TF2.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"

#include "Math/SpecFuncMathCore.h" 
#include "Math/SpecFuncMathMore.h" 

#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"

#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
using namespace std;
const double PI=3.1415926;


// Read data
AmplitudeLib N("/Users/suyoupeng/Downloads/rcbkdipole-master/data/proton/mve.dat");

int main(){
    return 0;
}
