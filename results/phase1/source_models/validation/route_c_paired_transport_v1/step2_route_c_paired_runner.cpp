#include <G4Box.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4Event.hh>
#include <G4Gamma.hh>
#include <G4LogicalVolume.hh>
#include <G4MTRunManager.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4ParticleGun.hh>
#include <G4Run.hh>
#include <G4RunManagerFactory.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4UImanager.hh>
#include <G4VModularPhysicsList.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VUserActionInitialization.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4Version.hh>
#include <G4ios.hh>
#include <Randomize.hh>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
constexpr int kRequiredThreads = 16;
constexpr long long kMinimumHistories = 100000000LL;

static std::string trim(std::string s) {
    auto ns=[](unsigned char c){return !std::isspace(c);};
    s.erase(s.begin(),std::find_if(s.begin(),s.end(),ns));
    s.erase(std::find_if(s.rbegin(),s.rend(),ns).base(),s.end());
    return s;
}
static std::vector<std::string> split(const std::string& s,char delim=',') {
    std::vector<std::string> out; std::stringstream ss(s); std::string item;
    while(std::getline(ss,item,delim)) out.push_back(trim(item));
    return out;
}
static std::map<std::string,std::string> readKv(const fs::path& p) {
    std::ifstream in(p); if(!in) throw std::runtime_error("Cannot open config: "+p.string());
    std::map<std::string,std::string> kv; std::string line;
    while(std::getline(in,line)) { line=trim(line); if(line.empty()||line[0]=='#') continue; auto pos=line.find('='); if(pos!=std::string::npos) kv[trim(line.substr(0,pos))]=trim(line.substr(pos+1)); }
    return kv;
}
static std::string req(const std::map<std::string,std::string>& kv,const std::string& k){auto it=kv.find(k);if(it==kv.end())throw std::runtime_error("Missing config key: "+k);return it->second;}
static double getd(const std::map<std::string,std::string>& kv,const std::string& k,double d=0){auto it=kv.find(k);return it==kv.end()?d:std::stod(it->second);}
static long long geti64(const std::map<std::string,std::string>& kv,const std::string& k,long long d=0){auto it=kv.find(k);return it==kv.end()?d:std::stoll(it->second);}

struct CsvTable {std::vector<std::string> header;std::vector<std::vector<std::string>> rows;std::map<std::string,size_t> idx;};
static CsvTable readCsv(const fs::path& p){
    std::ifstream in(p); if(!in) throw std::runtime_error("Cannot open CSV: "+p.string());
    CsvTable t;std::string line;if(!std::getline(in,line))throw std::runtime_error("Empty CSV: "+p.string());t.header=split(line);for(size_t i=0;i<t.header.size();++i)t.idx[t.header[i]]=i;
    while(std::getline(in,line)){if(trim(line).empty())continue;auto r=split(line);if(r.size()<t.header.size())r.resize(t.header.size());t.rows.push_back(std::move(r));}return t;
}
static std::string csvGet(const CsvTable& t,const std::vector<std::string>& r,const std::string& c){auto it=t.idx.find(c);if(it==t.idx.end())throw std::runtime_error("CSV missing column: "+c);return r.at(it->second);}

struct SourceBin{double low=0,high=0,prob=0;};
struct SourceSpectrum{
    std::vector<SourceBin> bins;std::vector<double> cdf;
    static SourceSpectrum load(const fs::path& p){
        CsvTable t=readCsv(p);SourceSpectrum s;double cum=0;
        for(const auto& r:t.rows){SourceBin b;b.low=std::stod(csvGet(t,r,"energy_low_MeV"));b.high=std::stod(csvGet(t,r,"energy_high_MeV"));
            if(t.idx.count("sampling_probability")) b.prob=std::stod(csvGet(t,r,"sampling_probability")); else b.prob=std::stod(csvGet(t,r,"probability_mass_bin"));
            if(!(b.high>b.low)||b.prob<0)throw std::runtime_error("Invalid source bin");cum+=b.prob;s.bins.push_back(b);s.cdf.push_back(cum);}
        if(s.bins.empty()||!(cum>0))throw std::runtime_error("Invalid source spectrum");for(auto& x:s.cdf)x/=cum;return s;
    }
    double sample() const {double u=G4UniformRand();auto it=std::lower_bound(cdf.begin(),cdf.end(),u);size_t i=std::min<size_t>(std::distance(cdf.begin(),it),bins.size()-1);const auto& b=bins[i];return std::max(1e-11,b.low+(b.high-b.low)*G4UniformRand());}
};

struct DepthBin{double center=0,low=0,high=0;};
static std::vector<DepthBin> loadDepthBins(const fs::path& p){
    CsvTable t=readCsv(p);std::vector<double> c;for(const auto& r:t.rows)c.push_back(std::stod(csvGet(t,r,"depth_cm")));std::sort(c.begin(),c.end());c.erase(std::unique(c.begin(),c.end()),c.end());
    if(c.size()<2||c.front()<0)throw std::runtime_error("PDD reference requires >=2 nonnegative depths");std::vector<double> e(c.size()+1);e[0]=std::max(0.0,c[0]-0.5*(c[1]-c[0]));
    for(size_t i=1;i<c.size();++i)e[i]=0.5*(c[i-1]+c[i]);e.back()=c.back()+0.5*(c.back()-c[c.size()-2]);std::vector<DepthBin> out;out.reserve(c.size());
    for(size_t i=0;i<c.size();++i){if(!(e[i+1]>e[i]))throw std::runtime_error("Invalid PDD depth bin width");out.push_back({c[i],e[i],e[i+1]});}return out;
}

struct Config{
    std::string runId,sourceModelId,sourcePlaneMode,phaseSpaceModel;
    double ssdCm=100,fieldXCm=40,fieldYCm=40,waterHalfXYCm=40,waterDepthCm=50,sourcePlaneOffsetCm=0.001;
    double pddScoreHalfWidthCm=1.0,profileDepthCm=10.0,profileDepthHalfWidthCm=0.05,profileYHalfWidthCm=1.0,profileBinWidthCm=0.1,profileHalfExtentCm=30.0;
    long long histories=0,seed=1;fs::path sourceCsv,pddReferenceCsv,pddOutputCsv,profileOutputCsv;
    static Config load(const fs::path& p){auto kv=readKv(p);Config c;c.runId=req(kv,"run_id");c.sourceModelId=req(kv,"source_model_id");c.ssdCm=getd(kv,"ssd_cm",100);c.fieldXCm=getd(kv,"field_x_cm",40);c.fieldYCm=getd(kv,"field_y_cm",40);c.waterHalfXYCm=getd(kv,"water_half_xy_cm",40);c.waterDepthCm=getd(kv,"water_depth_cm",50);c.sourcePlaneOffsetCm=getd(kv,"source_plane_offset_cm",0.001);c.pddScoreHalfWidthCm=getd(kv,"pdd_score_half_width_cm",1);c.profileDepthCm=getd(kv,"profile_depth_cm",10);c.profileDepthHalfWidthCm=getd(kv,"profile_depth_half_width_cm",0.05);c.profileYHalfWidthCm=getd(kv,"profile_y_half_width_cm",1);c.profileBinWidthCm=getd(kv,"profile_bin_width_cm",0.1);c.profileHalfExtentCm=getd(kv,"profile_half_extent_cm",30);c.histories=geti64(kv,"histories");c.seed=geti64(kv,"random_seed",1);c.sourceCsv=req(kv,"source_csv");c.pddReferenceCsv=req(kv,"pdd_reference_csv");c.pddOutputCsv=req(kv,"pdd_output_csv");c.profileOutputCsv=req(kv,"profile_output_csv");c.sourcePlaneMode=req(kv,"source_plane_mode");c.phaseSpaceModel=req(kv,"phase_space_model");
        if(c.histories<kMinimumHistories)throw std::runtime_error("Route-C paired run requires >=100,000,000 histories");if(!(c.ssdCm>0&&c.fieldXCm>0&&c.fieldYCm>0&&c.waterDepthCm>0&&c.sourcePlaneOffsetCm>0&&c.pddScoreHalfWidthCm>0&&c.profileDepthHalfWidthCm>0&&c.profileYHalfWidthCm>0&&c.profileBinWidthCm>0&&c.profileHalfExtentCm>0))throw std::runtime_error("Invalid Route-C geometry");if(c.sourcePlaneMode!="water_surface_factorized")throw std::runtime_error("Invalid source_plane_mode");if(c.phaseSpaceModel!="aggregate_energy_uniform_xy_virtual_source_divergence")throw std::runtime_error("Invalid phase_space_model");return c;}
};

class Detector final:public G4VUserDetectorConstruction{Config c_;public:explicit Detector(Config c):c_(std::move(c)){}G4VPhysicalVolume* Construct() override{auto*n=G4NistManager::Instance();auto*air=n->FindOrBuildMaterial("G4_AIR");auto*water=n->FindOrBuildMaterial("G4_WATER");double hxy=std::max(100.0,c_.waterHalfXYCm+10.0);double hz=std::max(30.0,c_.waterDepthCm+20.0);auto*ws=new G4Box("routec_world",hxy*cm,hxy*cm,hz*cm);auto*wl=new G4LogicalVolume(ws,air,"routec_world");auto*wp=new G4PVPlacement(nullptr,G4ThreeVector(),wl,"routec_world",nullptr,false,0);auto*ps=new G4Box("routec_water",c_.waterHalfXYCm*cm,c_.waterHalfXYCm*cm,0.5*c_.waterDepthCm*cm);auto*pl=new G4LogicalVolume(ps,water,"routec_water");new G4PVPlacement(nullptr,G4ThreeVector(0,0,0.5*c_.waterDepthCm*cm),pl,"routec_water",wl,false,0);return wp;}};

struct Accum{double sum=0,sumsq=0;void add(double x){sum+=x;sumsq+=x*x;}};
struct Shared{
    std::mutex mu;std::vector<Accum> pdd,profile;long long events=0;
    Shared(size_t np,size_t nx):pdd(np),profile(nx){}
    void merge(const std::vector<Accum>& a,const std::vector<Accum>& b,long long n){std::lock_guard<std::mutex>lk(mu);for(size_t i=0;i<pdd.size();++i){pdd[i].sum+=a[i].sum;pdd[i].sumsq+=a[i].sumsq;}for(size_t i=0;i<profile.size();++i){profile[i].sum+=b[i].sum;profile[i].sumsq+=b[i].sumsq;}events+=n;}
};

class RunAction;
class EventAction final:public G4UserEventAction{
    RunAction* run_;std::vector<double> pdd_,profile_;std::vector<size_t> tp_,tx_;
public:EventAction(RunAction*r,size_t np,size_t nx):run_(r),pdd_(np,0),profile_(nx,0){tp_.reserve(64);tx_.reserve(8);}void BeginOfEventAction(const G4Event*) override{for(auto i:tp_)pdd_[i]=0;for(auto i:tx_)profile_[i]=0;tp_.clear();tx_.clear();}void EndOfEventAction(const G4Event*) override;void addPdd(size_t i,double v){if(pdd_[i]==0)tp_.push_back(i);pdd_[i]+=v;}void addProfile(size_t i,double v){if(profile_[i]==0)tx_.push_back(i);profile_[i]+=v;}const std::vector<double>&pdd()const{return pdd_;}const std::vector<double>&profile()const{return profile_;}const std::vector<size_t>&tp()const{return tp_;}const std::vector<size_t>&tx()const{return tx_;}};
class RunAction final:public G4UserRunAction{
    Shared& s_;std::vector<Accum> pdd_,profile_;long long n_=0;
public:explicit RunAction(Shared&s):s_(s),pdd_(s.pdd.size()),profile_(s.profile.size()){}void addEvent(const std::vector<double>&p,const std::vector<size_t>&tp,const std::vector<double>&x,const std::vector<size_t>&tx){for(auto i:tp)pdd_[i].add(p[i]);for(auto i:tx)profile_[i].add(x[i]);++n_;}void EndOfRunAction(const G4Run*) override{s_.merge(pdd_,profile_,n_);}};
void EventAction::EndOfEventAction(const G4Event*){run_->addEvent(pdd_,tp_,profile_,tx_);}

class Primary final:public G4VUserPrimaryGeneratorAction{Config c_;SourceSpectrum src_;G4ParticleGun gun_{1};public:Primary(Config c,SourceSpectrum s):c_(std::move(c)),src_(std::move(s)){gun_.SetParticleDefinition(G4Gamma::Gamma());}void GeneratePrimaries(G4Event*e) override{gun_.SetParticleEnergy(src_.sample()*MeV);double x=(G4UniformRand()-0.5)*c_.fieldXCm;double y=(G4UniformRand()-0.5)*c_.fieldYCm;gun_.SetParticlePosition(G4ThreeVector(x*cm,y*cm,-c_.sourcePlaneOffsetCm*cm));G4ThreeVector d(x*cm,y*cm,c_.ssdCm*cm);gun_.SetParticleMomentumDirection(d.unit());gun_.GeneratePrimaryVertex(e);}};

class StepAction final:public G4UserSteppingAction{
    Config c_;const std::vector<DepthBin>& db_;std::vector<double> de_;EventAction& ev_;double profileEdgeMin_;size_t profileN_;
public:StepAction(Config c,const std::vector<DepthBin>&db,EventAction&e,size_t nx):c_(std::move(c)),db_(db),ev_(e),profileN_(nx){de_.reserve(db_.size()+1);de_.push_back(db_.front().low);for(const auto&b:db_)de_.push_back(b.high);profileEdgeMin_=-c_.profileHalfExtentCm-0.5*c_.profileBinWidthCm;}
    void UserSteppingAction(const G4Step*st) override{double edep=st->GetTotalEnergyDeposit()/MeV;if(!(edep>0))return;auto*pre=st->GetPreStepPoint();auto touch=pre->GetTouchableHandle();auto*vol=touch->GetVolume();if(!vol||vol->GetLogicalVolume()->GetName()!="routec_water")return;auto pos=0.5*(pre->GetPosition()+st->GetPostStepPoint()->GetPosition());double x=pos.x()/cm,y=pos.y()/cm,z=pos.z()/cm;if(z<0)return;
        if(std::abs(x)<=c_.pddScoreHalfWidthCm&&std::abs(y)<=c_.pddScoreHalfWidthCm&&z>=de_.front()&&z<=de_.back()){auto it=std::upper_bound(de_.begin(),de_.end(),z);size_t i=(it==de_.begin()?0:static_cast<size_t>(std::distance(de_.begin(),it)-1));if(i>=db_.size())i=db_.size()-1;double dz=db_[i].high-db_[i].low;double v=4.0*c_.pddScoreHalfWidthCm*c_.pddScoreHalfWidthCm*dz;if(v>0)ev_.addPdd(i,edep/v);}
        if(std::abs(y)<=c_.profileYHalfWidthCm&&std::abs(z-c_.profileDepthCm)<=c_.profileDepthHalfWidthCm){double u=(x-profileEdgeMin_)/c_.profileBinWidthCm;long long ii=static_cast<long long>(std::floor(u));if(ii>=0&&static_cast<size_t>(ii)<profileN_){double v=c_.profileBinWidthCm*(2*c_.profileYHalfWidthCm)*(2*c_.profileDepthHalfWidthCm);if(v>0)ev_.addProfile(static_cast<size_t>(ii),edep/v);}}
    }
};

class Actions final:public G4VUserActionInitialization{Config c_;SourceSpectrum s_;std::vector<DepthBin>d_;Shared& sh_;size_t nx_;public:Actions(Config c,SourceSpectrum s,std::vector<DepthBin>d,Shared&sh,size_t nx):c_(std::move(c)),s_(std::move(s)),d_(std::move(d)),sh_(sh),nx_(nx){}void BuildForMaster()const override{SetUserAction(new RunAction(sh_));}void Build()const override{auto*r=new RunAction(sh_);auto*e=new EventAction(r,sh_.pdd.size(),sh_.profile.size());SetUserAction(new Primary(c_,s_));SetUserAction(r);SetUserAction(e);SetUserAction(new StepAction(c_,d_,*e,nx_));}};
class Physics final:public G4VModularPhysicsList{public:Physics(){SetVerboseLevel(1);RegisterPhysics(new G4EmStandardPhysics_option4());}};

static std::pair<double,double> meanSem(const Accum&a,long long n){if(n<=0)return{0,0};double N=static_cast<double>(n),m=a.sum/N;if(n<2)return{m,0};double ss=std::max(0.0,a.sumsq-N*m*m);double sample=ss/(N-1.0);return{m,std::sqrt(sample/N)};}
static void writeOutputs(const Config&c,const std::vector<DepthBin>&db,const Shared&s){
    if(s.events!=c.histories)throw std::runtime_error("Merged event count differs from requested histories");std::vector<double> pm(db.size()),ps(db.size());for(size_t i=0;i<db.size();++i){auto q=meanSem(s.pdd[i],c.histories);pm[i]=q.first;ps[i]=q.second;}auto im=std::distance(pm.begin(),std::max_element(pm.begin(),pm.end()));if(!(pm.at(im)>0))throw std::runtime_error("PDD has zero normalization dose");
    fs::create_directories(c.pddOutputCsv.parent_path());std::ofstream po(c.pddOutputCsv);if(!po)throw std::runtime_error("Cannot open PDD output");po<<std::setprecision(17);po<<"run_id,source_model_id,depth_cm,mc_pdd_percent,mc_raw_mean_relative_dose,mc_raw_sem_relative_dose,normalization_depth_cm,histories,transport_threads,pdd_score_half_width_cm,ssd_cm,field_size_x_cm,field_size_y_cm,source_plane_mode,source_plane_offset_cm,phase_space_model,phase_space_complete,mc_uncertainty_note\n";for(size_t i=0;i<db.size();++i){double p=100.0*pm[i]/pm[im];po<<c.runId<<","<<c.sourceModelId<<","<<db[i].center<<","<<p<<","<<pm[i]<<","<<ps[i]<<","<<db[im].center<<","<<c.histories<<","<<kRequiredThreads<<","<<c.pddScoreHalfWidthCm<<","<<c.ssdCm<<","<<c.fieldXCm<<","<<c.fieldYCm<<","<<c.sourcePlaneMode<<","<<c.sourcePlaneOffsetCm<<","<<c.phaseSpaceModel<<",0,raw_event_sem_only_normalized_ratio_uncertainty_not_gating\n";}po.close();if(!po)throw std::runtime_error("Failed writing PDD output");
    const size_t nx=s.profile.size();const size_t cax=(nx-1)/2;std::vector<double> xm(nx),xs(nx);for(size_t i=0;i<nx;++i){auto q=meanSem(s.profile[i],c.histories);xm[i]=q.first;xs[i]=q.second;}if(!(xm.at(cax)>0))throw std::runtime_error("Profile CAX dose is zero");fs::create_directories(c.profileOutputCsv.parent_path());std::ofstream xo(c.profileOutputCsv);if(!xo)throw std::runtime_error("Cannot open profile output");xo<<std::setprecision(17);xo<<"run_id,source_model_id,distance_cm,mc_relative_dose_percent,mc_raw_mean_relative_dose,mc_raw_sem_relative_dose,normalization_distance_cm,profile_depth_cm,profile_depth_half_width_cm,profile_y_half_width_cm,profile_bin_width_cm,histories,transport_threads,ssd_cm,field_size_x_cm,field_size_y_cm,source_plane_mode,source_plane_offset_cm,phase_space_model,phase_space_complete,mc_uncertainty_note\n";for(size_t i=0;i<nx;++i){double x=-c.profileHalfExtentCm+static_cast<double>(i)*c.profileBinWidthCm;double p=100.0*xm[i]/xm[cax];xo<<c.runId<<","<<c.sourceModelId<<","<<x<<","<<p<<","<<xm[i]<<","<<xs[i]<<",0,"<<c.profileDepthCm<<","<<c.profileDepthHalfWidthCm<<","<<c.profileYHalfWidthCm<<","<<c.profileBinWidthCm<<","<<c.histories<<","<<kRequiredThreads<<","<<c.ssdCm<<","<<c.fieldXCm<<","<<c.fieldYCm<<","<<c.sourcePlaneMode<<","<<c.sourcePlaneOffsetCm<<","<<c.phaseSpaceModel<<",0,raw_event_sem_only_normalized_ratio_uncertainty_not_gating\n";}xo.close();if(!xo)throw std::runtime_error("Failed writing profile output");
}

static int runPaired(const fs::path&cfgPath){
#ifndef G4MULTITHREADED
    throw std::runtime_error("Route-C runner requires multithreaded Geant4");
#else
    Config c=Config::load(cfgPath);SourceSpectrum src=SourceSpectrum::load(c.sourceCsv);auto db=loadDepthBins(c.pddReferenceCsv);const size_t nx=static_cast<size_t>(std::llround(2.0*c.profileHalfExtentCm/c.profileBinWidthCm))+1;if(nx<3||nx%2==0)throw std::runtime_error("Profile grid must contain odd number of bins including CAX");Shared sh(db.size(),nx);G4Random::setTheSeed(static_cast<long>(c.seed));auto*rm=G4RunManagerFactory::CreateRunManager(G4RunManagerType::MTOnly);auto*mt=dynamic_cast<G4MTRunManager*>(rm);if(!mt)throw std::runtime_error("Geant4 MT run manager unavailable");mt->SetNumberOfThreads(kRequiredThreads);rm->SetUserInitialization(new Detector(c));rm->SetUserInitialization(new Physics());rm->SetUserInitialization(new Actions(c,src,db,sh,nx));auto*ui=G4UImanager::GetUIpointer();ui->ApplyCommand("/control/verbose 1");ui->ApplyCommand("/run/verbose 1");ui->ApplyCommand("/event/verbose 0");ui->ApplyCommand("/tracking/verbose 0");G4cout<<"STEP2_ROUTE_C_PAIRED_RUNNER=1"<<G4endl;G4cout<<"WORKER_THREADS="<<kRequiredThreads<<G4endl;G4cout<<"HISTORIES="<<c.histories<<G4endl;G4cout<<"PROFILE_BIN_WIDTH_CM="<<c.profileBinWidthCm<<G4endl;G4cout<<"PROFILE_DEPTH_HALF_WIDTH_CM="<<c.profileDepthHalfWidthCm<<G4endl;G4cout<<"PROFILE_Y_HALF_WIDTH_CM="<<c.profileYHalfWidthCm<<G4endl;rm->Initialize();rm->BeamOn(static_cast<G4int>(c.histories));writeOutputs(c,db,sh);delete rm;return 0;
#endif
}

int main(int argc,char**argv){try{if(argc<2){std::cerr<<"usage: step2_route_c_runner info | paired <config.ini>\n";return 2;}std::string m=argv[1];if(m=="info"){std::cout<<"IMPLEMENTATION_LANGUAGE=C++17\nREQUIRED_THREADS="<<kRequiredThreads<<"\nMINIMUM_HISTORIES="<<kMinimumHistories<<"\nGEANT4_VERSION="<<G4Version<<"\n";
#ifdef G4MULTITHREADED
std::cout<<"MULTITHREADED=1\n";
#else
std::cout<<"MULTITHREADED=0\n";
#endif
return 0;}if(m=="paired"&&argc==3)return runPaired(argv[2]);throw std::runtime_error("Invalid arguments");}catch(const std::exception&e){std::cerr<<"FATAL: "<<e.what()<<"\n";return 1;}}
