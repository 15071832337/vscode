#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>
#include <TMath.h>

void ReadData() {
    // 打开文件
    TFile* file = TFile::Open("/home/gaoyh/Downloads/HFJetTreeProduction/output/hf_jetsTree_0.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }
    
    // 获取树 - 使用正确的名称
    TTree* tree = (TTree*)file->Get("hfTree");
    if (!tree) {
        std::cerr << "Error: Could not find tree 'hfTree'" << std::endl;
        file->Close();
        return;
    }
    // Print all branch names
    TObjArray* branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* branch = (TBranch*)branches->At(i);
        std::cout << "Branch name: " << branch->GetName() << std::endl;
    }
    // 声明变量 - 使用指针
    int multiplicity = 0;
    int fEventId = 0;
    std::vector<int>* particlesid = nullptr;
    std::vector<int>* hadronPDG = nullptr;
    std::vector<double>* hadronPt = nullptr;
    std::vector<double>* fHadronEta = nullptr;
    std::vector<double>* fHadronPhi = nullptr;
    std::vector<double>* fHadronPx = nullptr;
    std::vector<double>* fHadronPy = nullptr;
    std::vector<double>* fHadronPz = nullptr;
    std::vector<double>* fHadronE = nullptr;
    std::vector<bool>* fHadronIsCharged = nullptr;
    // ... 添加其他需要的向量 ...
    
    // 设置分支地址 - 使用指针的地址
    tree->SetBranchAddress("eventId", &fEventId);
    tree->SetBranchAddress("multiplicity", &multiplicity);
    tree->SetBranchAddress("particlesid", &particlesid);
    tree->SetBranchAddress("hadronPDG", &hadronPDG);
    tree->SetBranchAddress("hadronPt", &hadronPt);

    tree->SetBranchAddress("hadronEta", &fHadronEta);
    tree->SetBranchAddress("hadronPhi", &fHadronPhi);
    tree->SetBranchAddress("hadronPx", &fHadronPx);
    tree->SetBranchAddress("hadronPy", &fHadronPy);
    tree->SetBranchAddress("hadronPz", &fHadronPz);
    tree->SetBranchAddress("hadronE", &fHadronE);
    tree->SetBranchAddress("hadronIsCharged", &fHadronIsCharged);

    // ... 为其他分支设置地址 ...
    
    // 创建直方图
    TH1F* hhadronPt = new TH1F("hhadronPt", "Particle pT Distribution;pT (GeV/c);Counts", 100, 0, 10);
    TH2D* hhadroneta_vs_hadronphi = new TH2D("hhadroneta_vs_hadronphi", "hadroneta vs hadronphi;#eta;#phi", 100, -5, 5, 100, 0, 2 * TMath::Pi());
    
    // 遍历事件
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Total events: " << nEntries << std::endl;
    
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        // 检查向量是否有效
        if (!particlesid || !hadronPDG || !hadronPt) {
            std::cerr << "Error: Null vector pointer in event " << i << std::endl;
            continue;
        }
        
        // 检查大小一致性
        if (multiplicity != (int)particlesid->size() || 
            multiplicity != (int)hadronPDG->size() ||
            multiplicity != (int)hadronPt->size() ||
            multiplicity != (int)fHadronEta->size() ||
            multiplicity != (int)fHadronPhi->size()) {
            std::cerr << "Warning: Size mismatch in event " << fEventId
                      << " - multiplicity: " << multiplicity
                      << ", particlesid: " << particlesid->size()
                      << ", fHadronEta: " << fHadronEta->size()
                      << ", fHadronPhi: " << fHadronPhi->size()
                      << ", hadronPDG: " << hadronPDG->size()
                      << ", hadronPt: " << hadronPt->size() << std::endl;
        }
        
        std::cout << "\nEvent ID: " << fEventId
                  << ", Particles: " << multiplicity << std::endl;
        
        // 遍历当前事件的所有粒子
        for (int j = 0; j < multiplicity; j++) {
            // 检查索引是否在范围内
            if (j >= particlesid->size() || j >= hadronPDG->size() || j >= hadronPt->size()) {
                std::cerr << "Error: Index " << j << " out of range in event " << fEventId << std::endl;
                break;
            }
            
            int pid = particlesid->at(j);
            int pdg = hadronPDG->at(j);
            double pt = hadronPt->at(j);
            double eta = fHadronEta->at(j);
            double phi = fHadronPhi->at(j);
            
            std::cout << "  Particle " << j << ": "
                      << "ID = " << pid << ", "
                      << "PDG = " << pdg << ", "
                      << "pT = " << pt << ", "
                      << "#eta = " << eta << ", "
                      << "#phi = " << phi << std::endl;

            
            // 填充直方图 - 添加基本检查
            if (pt > 0 && pt < 1000) {  // 合理的pT范围
                hhadronPt->Fill(pt);
            } else {
                std::cerr << "Skipping invalid pT value: " << pt << " in event " << fEventId << std::endl;
            }
            // 填充eta-phi图
            /*if (eta >= -5 && eta <= 5 && phi >= 0 && phi <= 2 * TMath::Pi()) {
                hhadroneta_vs_hadronphi->Fill(eta, phi);
            } else {
                std::cerr << "Skipping invalid eta or phi value: eta = " << eta << ", phi = " << phi << " in event " << fEventId << std::endl;
            }*/
            hhadroneta_vs_hadronphi->Fill(eta, phi);
        }
    }
    
    // 绘制直方图
    TCanvas* c1 = new TCanvas("c1", "pT Distribution", 800, 600);
    hhadronPt->Draw();
    c1->Update();
    
    // 绘制eta-phi图
    TCanvas* c2 = new TCanvas("c2", "eta vs phi", 800, 600);
    hhadroneta_vs_hadronphi->Draw("colz");  // "colz"选项显示颜色表示密度
    c2->Update();  
    c2->Draw();

    // 保持窗口打开（仅用于交互式会话）
    // c1->WaitPrimitive();
    
    // 清理
    // file->Close();  // 保持文件打开以便访问直方图
    // 注意：不要删除ROOT管理的向量指针
}