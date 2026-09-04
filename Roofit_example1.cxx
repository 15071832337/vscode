//ROOT provides with the RooFit library a toolkit for modeling the expected distribution
// of events in a physics analysis. Models can be used to perform unbinned maximum likelihood fits, 
//create plots, and generate “toy Monte Carlo” samples for various studies.

/*
RooFit的核心功能是支持对“事件数据”分布进行建模，其中每个事件是时间上的一个离散发生，并且具有一个或多个与之相关的测量观测量。
这类实验产生的数据集遵循泊松（或二项式）统计规律。对此类分布进行自然建模的语言是概率密度函数（PDF）F(x;p)，它描述了观测量x的分布概率密度，
表示为参数p的函数。
概率密度函数的定义性质——在所有观测量上的单位归一化以及正定性——也为结构化建模语言的设计提供了重要优势：PDF可以很容易地相加，
且分数系数具有直观的物理含义。它们允许用低维基本模块构建高维PDF，并采用直观的语言来引入和描述观测量之间的相关性。
同时，它们也使得玩具蒙特卡洛抽样技术的通用实现成为可能，并且当然也是使用（非分箱）最大似然参数估计方法的先决条件。

RooFit在其数学数据模型组件到C++对象的映射中引入了一种细粒度结构：其目标并非创建一个描述数据模型的单一整体实体，而是将每个数学符号都表示为独立的对象。这种设计理念的一个特点是，所有RooFit模型始终由多个对象构成。
数学概念 RooFit类

变量	RooRealVar
函数	RooAbsReal
概率密度函数（PDF）	RooAbsPdf
积分	RooRealIntegral
空间点	RooArgSet
空间点列表	RooAbsData


*/

// 或

void Roofit_example1(){
    using namespace RooFit;

    TFile *outFile = new TFile("argus_model_plot.root", "RECREATE");
    TCanvas *c = new TCanvas("c1", "Signal + Background", 1000, 700);
  
    // ===== 1. 观测量 =====
    RooRealVar mes("mes", "m_{ES} (GeV)", 5.20, 5.30);

    // ===== 2. 信号参数 =====
    RooRealVar sigmean("sigmean", "B^{#pm} mass", 5.28, 5.20, 5.30);
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
    //(nsig*signal+nbkg*linearBg)/(nsig+nbkg)

    // ===== 5. 生成数据 =====
    std::unique_ptr<RooDataSet> data{model.generate(mes, 20000)};

    // ===== 6. 拟合 =====
    RooFitResult* result = model.fitTo(*data, Save(true));

    // ===== 7. 绘图 =====
    RooPlot *mesframe = mes.frame();
    data->plotOn(mesframe);
    model.plotOn(mesframe, LineColor(kBlack), LineWidth(10));
    model.plotOn(mesframe, Components(linearBg), LineStyle(ELineStyle::kDashed), 
                 LineColor(kBlue), LineWidth(5));
    model.plotOn(mesframe, Components(signal), LineColor(kRed), 
                 LineStyle(ELineStyle::kDotted), LineWidth(2));

    mesframe->SetTitle("Signal + Linear Background;m_{ES} (GeV);Events / (0.001)");

    // ====================================================
    // ===== 8. 在图上显示拟合参数（先画图，再画文本） =====
    // ====================================================
    
    // ★ 先画 RooPlot
    mesframe->Draw();

    // ★ 然后在上面叠加 TPaveText
    TPaveText *pt = new TPaveText(0.35, 0.30, 0.65, 0.85, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.030);
    pt->SetTextFont(42);

    pt->AddText("--- Fit Results ---");
    pt->AddText(" ");
    pt->AddText(Form("N_{sig}  = %.1f #pm %.1f", nsig.getVal(), nsig.getError()));
    pt->AddText(Form("N_{bkg}  = %.1f #pm %.1f", nbkg.getVal(), nbkg.getError()));
    pt->AddText(" ");
    pt->AddText(Form("#mu      = %.5f #pm %.5f", sigmean.getVal(), sigmean.getError()));
    pt->AddText(Form("#sigma   = %.5f #pm %.5f", sigwidth.getVal(), sigwidth.getError()));
    pt->AddText(" ");
    pt->AddText(Form("c0       = %.3f #pm %.3f", c0.getVal(), c0.getError()));
    pt->AddText(Form("c1       = %.3f #pm %.3f", c1.getVal(), c1.getError()));

    if (result) {
        pt->AddText(" ");
        pt->AddText(Form("Status   = %d", result->status()));
        pt->AddText(Form("Cov Qual = %d", result->covQual()));
    }

    pt->Draw();
    c->Update();

    // ===== 9. 保存 =====
    outFile->cd();
    mesframe->Write("mes_frame");
    c->Write("c1");
    c->SaveAs("fit_result.png");
    
    // 打印到终端
    std::cout << "\n========== 拟合结果 ==========" << std::endl;
    result->Print("v");
    
    delete c;
    outFile->Close();
}














