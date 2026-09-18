//ROOT provides with the RooFit library a toolkit for modeling the expected distribution
// of events in a physics analysis. Models can be used to perform unbinned maximum likelihood fits, 
//create plots, and generate “toy Monte Carlo” samples for various studies.

/*
bin counting 技术

*/

// 或

Bool_t rejectCenter = kTRUE;
Double_t g_xMin = 0;
Double_t g_xMax = 0;

// ===== 自定义拟合函数 =====
Double_t myLinearFit(Double_t *x, Double_t *par)
{
    // 如果启用拒绝，且 x 在中心区域内，则拒绝这个点
    
    if (rejectCenter && x[0] > g_xMin && x[0] < g_xMax) {
        TF1::RejectPoint();  // 关键：调用静态函数拒绝该点
        return 0;
    }

    // 线性函数：par[0] + par[1]*x
    return par[0] + par[1] * x[0];
}

void Roofit_example2(){
    using namespace RooFit;

    TFile *outFile = new TFile("bincounting.root", "RECREATE");
    TCanvas *c = new TCanvas("c1", "Signal + Background", 2500, 1000);
   
    // ===== 1. 观测量 =====
    RooRealVar mes("mes", "m_{ES} (GeV)", 5.20, 5.30);

    // ===== 2. 信号参数 =====
    RooRealVar sigmean("sigmean", "B^{#pm} mass", 5.25, 5.20, 5.30);
    RooRealVar sigwidth("sigwidth", "B^{#pm} width", 0.0027, 0.001, 1.0);
    RooGaussian signal("signal", "signal PDF", mes, sigmean, sigwidth);

    // ===== 3. 背景参数 =====
    RooRealVar c0("c0", "intercept", 478.0, 0.0, 1000.0);
    RooRealVar c1("c1", "slope", -90.0, -200.0, -1.0);
    RooPolynomial linearBg("background", "Linear background", mes, RooArgList(c0, c1));

    // ===== 4. 复合模型 =====
    RooRealVar nsig("nsig", "#signal events", 500, 0., 20000);
    RooRealVar nbkg("nbkg", "#background events", 1500, 0., 20000);
    RooAddPdf model("model", "signal + background", {signal, linearBg}, {nsig, nbkg});

    // ===== 5. 生成数据 =====
    std::unique_ptr<RooDataSet> data{model.generate(mes, 20000)};

    // ===== 6. 拟合 =====
    RooFitResult* result = model.fitTo(*data, Save(true));

    // ===== 7. 绘图 =====
    RooPlot *mesframe = mes.frame();
    data->plotOn(mesframe);
    model.plotOn(mesframe, LineColor(kBlack), LineWidth(2));
    //model.plotOn(mesframe, Components(linearBg), LineStyle(ELineStyle::kDashed), LineColor(kBlue), LineWidth(2));
   // model.plotOn(mesframe, Components(signal), LineColor(kRed), LineStyle(ELineStyle::kDotted), LineWidth(2));


    // ========================================================
    // ===== 8. 在图上画出 ±6σ 的信号区域 =====
    // ========================================================

    RooHist *hist = dynamic_cast<RooHist*>(mesframe->getObject(0));
    int nBins = 60;  // 默认值
    
    nBins = hist->GetN();  // 获取实际分箱数
    std::cout << "RooPlot 分箱数: " << nBins << std::endl;
    
    
    mesframe->SetTitle("Signal + Linear Background;m_{ES} (GeV);Events / (0.001)");

    
      // 获取拟合后的均值和宽度
      double mean = sigmean.getVal();
      double sigma = sigwidth.getVal();
      
      // 计算 ±6σ 边界
      double xMin = mean - 6 * sigma;
      double xMax = mean + 6 * sigma;
  
      double xMin12_sigma = mean - 12* sigma;
      double xMax12_sigma = mean + 12* sigma;
      
      // 获取当前画布的坐标范围
      double yMin = mesframe->GetMinimum();
      double yMax = mesframe->GetMaximum();
  
  
      //bin counting method
      // 固定边界参数
      RooRealVar mu_minus_6sig("mu_minus_6sig", "mu-6sig", mean - 6*sigma);
      RooRealVar mu_plus_6sig("mu_plus_6sig", "mu+6sig", mean + 6*sigma);
      mu_minus_6sig.setConstant(true);
      mu_plus_6sig.setConstant(true);
  

      // ===== 使用 RooFit::Range() 指定多个区间 =====
  
  // ===== 共享参数的线性函数 =====
      RooRealVar a("a", "slope", 0.0, -1e-6, 1e-6);
      RooRealVar b("b", "intercept", 1.5e-8, 0.0, 1e-6);
      RooGenericPdf tailPDF(
      "tailPDF",
      "a * mes + b",
      RooArgList(mes, a, b)
      );
  
      std::unique_ptr<RooDataSet> leftData{
          static_cast<RooDataSet*>(data->reduce(RooFit::Range(xMin12_sigma, xMin)))
      };
      // 右尾数据
      std::unique_ptr<RooDataSet> rightData{
          static_cast<RooDataSet*>(data->reduce(RooFit::Range(xMax, xMax12_sigma)))
      };
      
      // 合并尾部数据
      RooDataSet tailData("tailData", "Tail Data", RooArgSet(mes));
  
      // 添加左尾数据
      for (int i = 0; i < leftData->numEntries(); i++) {
          const RooArgSet* row = leftData->get(i);
          mes.setVal(row->getRealValue("mes"));
          tailData.add(RooArgSet(mes));
      }
  
      // 添加右尾数据
      for (int i = 0; i < rightData->numEntries(); i++) {
          const RooArgSet* row = rightData->get(i);
          mes.setVal(row->getRealValue("mes"));
          tailData.add(RooArgSet(mes));
      }
      
      // 用共享参数的线性函数拟合合并后的数据
      RooFitResult* resultTail = tailPDF.fitTo(tailData, Save(true));


    // ===== 定义三个范围 =====
    std::string leftRange = Form("leftRange_%.4f_%.4f", xMin12_sigma, xMin);
    std::string rightRange = Form("rightRange_%.4f_%.4f", xMax, xMax12_sigma);
 

    // 在 mes 上定义范围
    mes.setRange(leftRange.c_str(), xMin12_sigma, xMin);
    mes.setRange(rightRange.c_str(), xMax, xMax12_sigma);



        // ============================================================
    // ★ 直接用 data 拟合，拒绝中心区域
    // ============================================================

    // 创建直方图
   
    g_xMin = xMin;
    g_xMax = xMax;
    rejectCenter = kTRUE;
      //xMin12_sigma，xMax12_sigma
    TH1D *hData = new TH1D("hData", "Full Data", nBins, 5.2,5.3);
    for (int i = 0; i < data->numEntries(); i++) {
        const RooArgSet* row = data->get(i);
        hData->Fill(row->getRealValue("mes"));
    }
    
    // ============================================================
    // ★ 创建 TF1 并传入自定义拟合函数
    // ============================================================
    TF1 *fitFunc = new TF1("fitFunc", myLinearFit, xMin12_sigma, xMax12_sigma, 2);
    fitFunc->SetParameters(1000, 0);  // 初始值：斜率=1，截距=0
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineStyle(ELineStyle::kDashed);

    // ★ 拟合（"Q" 表示安静模式）
    hData->Fit(fitFunc, "Q");

    double a_real = fitFunc->GetParameter(0);  // 截距
    double b_real = fitFunc->GetParameter(1);  // 斜率

    std::cout << "========== 拟合结果 ==========" << std::endl;
    std::cout << "斜率 b = " << b_real << std::endl;
    std::cout << "截距 a = " << a_real << std::endl;
    
    TF1 *fitFuncClone = (TF1*)fitFunc->Clone("fitFuncClone");
    
    // ★ 先画 RooPlot
    
    mesframe->Draw();

    fitFuncClone->Draw("same");
       

    TBox *box = new TBox(xMin, yMin, xMax, yMax);
    box->SetFillColor(kGreen);
    box->SetFillStyle(3001);  // 半透明样式
    box->SetLineColor(0);      // 无边框
    box->Draw();

    TBox *box_right = new TBox(xMax, yMin, xMax12_sigma, yMax);
    box_right->SetFillColor(kBlue);
    box_right->SetFillStyle(3001);  // 半透明样式
    box_right->SetLineColor(0);      // 无边框
    box_right->Draw();


    TBox *box_left = new TBox(xMin12_sigma, yMin, xMin, yMax);
    box_left->SetFillColor(kBlue);
    box_left->SetFillStyle(3001);  // 半透明样式
    box_left->SetLineColor(0);      // 无边框
    box_left->Draw();

    TPaveText *pt = new TPaveText(0.15, 0.40, 0.35, 0.85, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.030);
    pt->SetTextFont(42);
   

    pt->AddText("--- Fit Results ---");

    pt->AddText(Form("#mu      = %.5f #pm %.5f", sigmean.getVal(), sigmean.getError()));
    pt->AddText(Form("#sigma   = %.5f #pm %.5f", sigwidth.getVal(), sigwidth.getError()));
    pt->AddText(Form("c0       = %.3f #pm %.3f", c0.getVal(), c0.getError()));
    pt->AddText(Form("c1       = %.3f #pm %.3f", c1.getVal(), c1.getError()));
    pt->AddText(Form("a       = %.10f #pm %.15f", a.getVal(), a.getError()));
    pt->AddText(Form("b       = %.10f #pm %.15f", b.getVal(), b.getError()));

    pt->Draw();
    c->Update();    
    

    // ===== 9. 保存 =====
    outFile->cd();
    hData->Write("hdata");
    mesframe->Write("mes_frame");
    c->Write("c1");
    c->SaveAs("fit_result.png");
    
    // 打印到终端
    std::cout << "\n========== 拟合结果 ==========" << std::endl;
    result->Print("v");
    
    delete c;
    outFile->Close();


}














