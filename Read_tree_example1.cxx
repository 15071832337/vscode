/*
TTree的行为就像存储上的数据结构数组——除了一个条目（或数据库语言中的行）。该条目可以在内存中访问：
您可以加载任何树条目，最好是按顺序加载。您可以以变量的形式，為當前條目的列值提供您自己的存储空間。
在这种情况下，您必须告诉TTree这些变量的地址；要么通过调用TTree::SetBranchAddress()要么在创建分支进行写入时传递变量。
当“填充”（写入）TTree时，它将从这些变量中读取值；当读回aTTree条目时，它将从存储中读取的值写入您的变量中。
*/

/*
树由一个独立的列列表组成，称为分支。分支可以包含任何基本类型的值、ROOT类型系统已知的C++对象或这些对象的集合。
读取树时，您可以选择应该读取的分支子集。这允许您优化给定分析的读取吞吐量，并且是以列式格式存储数据的主要动机之一。
分支由TBranch及其派生类表示。虽然TBranch表示结构，但从TLeaf继承的对象可以访问实际数据。
最初，任何列数据都可以通过TLeaf访问；如今，一些TBranch衍生的类本身提供数据访问，如TBranchElement。
*/

void write_tree(){

    //编写TTree时，您首先要创建TFile（请参阅→ROOT文件。然后构建要存储在文件中的TTree；我们稍后将向树中添加分支。
        std::unique_ptr<TFile> myFile( TFile::Open("file.root", "RECREATE") ); 
        auto tree = std::make_unique<TTree>("tree", "The Tree Title");
    
    //如果您有int、float、bool或任何其他基本类型的变量类型，您可以从中创建一个分支（和叶子）。
    //对于基本数据类型，可以从变量中推断类型，叶子的名称将设置为分支的名称。
        float var;
        tree->Branch("branch0", &var);
    //"branch0"是这个树叶的名字
    
    //填满一棵树永久链接使用TTree:Fill()向树添加新条目（或“行”），并存储在创建分支期间提供的变量的当前值。
    //写树头永久链接,使用TTree::Write()将树头写入ROOT文件。早期条目的数据可能已经写入TTree::Fill()的一部分。
    
    for (int iEntry = 0; iEntry < 1000; ++iEntry) {
        var = 0.3 * iEntry;
        // Fill the current value of `var` into `branch0`
        tree->Fill();
     }
     //你还可以添加其他树叶
     double var1;
     tree->Branch("branch1", &var1);
    
     for (int iEntry = 0; iEntry < 10; ++iEntry) {
        var1 = iEntry*iEntry;
        // Fill the current value of `var` into `branch0`
        tree->Fill();
     }
    
     // Now write the header
     tree->Write();
    
    }
    
    void read_tree(){
        std::unique_ptr<TFile> myFile( TFile::Open("/Users/suyoupeng/Downloads/Roofit_study/file.root") );
        auto tree = myFile->Get<TTree>("tree");
    
        float variable;
        tree->SetBranchAddress("branch0", &variable);
    
        for (int iEntry = 0; iEntry<= 30; ++iEntry) {
            // Load the data for the given tree entry
            tree->GetEntry(iEntry);
    
            // Now, `variable` is set to the value of the branch
            // "branchName" in tree entry `iEntry`
            printf("%f\n", variable);
        }
        tree->Print();
        tree->Scan();
    
    }
    /*附加 TTree 作为 TChain 永久链接
    在高能物理中，你总是想要尽可能多的数据。但处理数 TB 级别的文件并不是件愉快的事。
    ROOT 允许你将数据分散存储到多个文件中，然后你可以将这些文件中各自的树部分作为一个大数来访问。
    这可以通过 TChain 来实现，TChain 继承自 TTree：它需要知道文件中树的名称（在添加文件时可以覆盖）以及文件名，
    然后它会表现得像一棵巨大且连续的树：
    
    TChain chain("CommonTreeName");
    if (chain.Add("data_*.root") != 12)
       std::cerr << "Expected to find 12 files!\n";
    // Use `chain` as if it was a `TTree`
    
    */
    
    
    void createMainTree() {
        // 创建输出文件
        TFile *f = new TFile("main.root", "RECREATE");
        TTree *tree = new TTree("tree", "Main tree with raw data");
    
        // ===== 定义分支 =====
        int event;
        float pt, eta, phi;
        int charge;
        tree->Branch("event", &event, "event/I");
        tree->Branch("pt", &pt, "pt/F");
        tree->Branch("eta", &eta, "eta/F");
        tree->Branch("phi", &phi, "phi/F");
        tree->Branch("charge", &charge, "charge/I");
    
        // ===== 填充数据 =====
        // 为了演示，手动填充 3 个事件
        // 实际使用时可以从数据文件读取或随机生成
        
        // 事件 1
        event = 1; pt = 45.0; eta = 0.5; phi = 2.1; charge = 1;
        tree->Fill();
        
        // 事件 2
        event = 2; pt = 78.0; eta = -0.3; phi = 0.8; charge = -1;
        tree->Fill();
        
        // 事件 3
        event = 3; pt = 23.0; eta = 1.2; phi = 3.4; charge = 1;
        tree->Fill();
        tree->Scan();
    
        // 保存并关闭
        tree->Write();
        f->Close();
        delete f;
        
        std::cout << "✅ main.root 创建完成，共 3 个事件" << std::endl;
    }
    
    void createFriendTree() {
        // 创建输出文件
        TFile *f = new TFile("friend.root", "RECREATE");
        TTree *tree = new TTree("friendTree", "Friend tree with corrected/calculated quantities");
    
        // ===== 定义分支 =====
        int event;
        float pt_corr, mass;
        int isGood;
        tree->Branch("event", &event, "event/I");
        tree->Branch("pt_corr", &pt_corr, "pt_corr/F");
        tree->Branch("mass", &mass, "mass/F");
        tree->Branch("isGood", &isGood, "isGood/I");
    
        // ===== 填充数据（事件顺序必须与 main.root 一致） =====
        // 事件 1
        event = 1; pt_corr = 45.2; mass = 5.2; isGood = 1;
        tree->Fill();
        
        // 事件 2
        event = 2; pt_corr = 77.5; mass = 4.8; isGood = 0;
        tree->Fill();
        
        // 事件 3
        event = 3; pt_corr = 23.1; mass = 5.5; isGood = 1;
        tree->Fill();
        tree->Scan();
        // 保存并关闭
        tree->Write();
        f->Close();
        delete f;
        
        std::cout << "✅ friend.root 创建完成，共 3 个事件" << std::endl;
    }
    
    void mergeAndView() {
        // ★ 使用 TChain
        TChain chain("tree");
        chain.Add("main.root");
        
        // ★ 添加友元树
        chain.AddFriend("friendTree", "friend.root");
        
        std::cout << "总事件数: " << chain.GetEntries() << std::endl;
        chain.Print("branches");
        chain.Show(1);
        
        int event, charge, isGood;
        float pt, eta, phi, pt_corr, mass;
        
        chain.SetBranchAddress("event", &event);
        chain.SetBranchAddress("pt", &pt);
        chain.SetBranchAddress("eta", &eta);
        chain.SetBranchAddress("phi", &phi);
        chain.SetBranchAddress("charge", &charge);
        
        chain.SetBranchAddress("friendTree.pt_corr", &pt_corr);
        chain.SetBranchAddress("friendTree.mass", &mass);
        chain.SetBranchAddress("friendTree.isGood", &isGood);
        
        
        std::cout << "\n========== 合并后的数据 ==========" << std::endl;
        std::cout << "event\tpt\teta\tphi\tcharge\tpt_corr\tmass\tisGood" << std::endl;
        std::cout << "--------------------------------------------------------" << std::endl;
        
        for (int i = 0; i < chain.GetEntries(); i++) {
            chain.GetEntry(i);
            printf("%d\t%.1f\t%.1f\t%.1f\t%d\t%.1f\t%.1f\t%d\n", 
                   event, pt, eta, phi, charge, pt_corr, mass, isGood);
        }
       
    }
    
    // ★ 运行所有步骤
    void runExample() {
        createMainTree();
        createFriendTree();
        mergeAndView();
    }



    void writeVariableParticleTreeRandom() {
        TFile *f = new TFile("particles_random.root", "RECREATE");
        TTree *tree = new TTree("tree", "Tree with variable number of particles");
    
        // ===== 定义分支变量 =====
        int eventNumber;
        std::vector<float> *pt = new std::vector<float>();
        std::vector<float> *eta = new std::vector<float>();
        std::vector<float> *phi = new std::vector<float>();
        std::vector<int> *charge = new std::vector<int>();
    
        tree->Branch("event", &eventNumber, "event/I");
        tree->Branch("pt", &pt);
        tree->Branch("eta", &eta);
        tree->Branch("phi", &phi);
        tree->Branch("charge", &charge);
    
        // ===== 随机数生成器 =====
        TRandom3 rnd(0);
    
        // ===== 生成 10 个事件 =====
        int nEvents = 10;
        
        for (int iEvent = 0; iEvent < nEvents; iEvent++) {
            eventNumber = iEvent + 1;  // 事件编号从 1 开始
            
            // ★ 随机决定这个事件有多少个粒子（1 到 10 个）
            int nParticles = rnd.Integer(10) + 1;
            
            // 清空 vector
            pt->clear();
            eta->clear();
            phi->clear();
            charge->clear();
            
            // ★ 填充粒子数据
            for (int iParticle = 0; iParticle < nParticles; iParticle++) {
                pt->push_back(rnd.Uniform(0, 100));           // pt: 0-100 GeV
                eta->push_back(rnd.Uniform(-3, 3));           // eta: -3 到 3
                phi->push_back(rnd.Uniform(-TMath::Pi(), TMath::Pi()));  // phi: -π 到 π
                charge->push_back(rnd.Uniform(0, 1) > 0.5 ? 1 : -1);    // charge: ±1
            }
            
            tree->Fill();
        }
    
        // ===== 打印树信息 =====
        tree->Print();
        tree->Scan();
    
        tree->Write();
        f->Close();
        delete f;
        
        std::cout << "✅ particles_random.root 创建完成，共 " << nEvents << " 个事件" << std::endl;
    }