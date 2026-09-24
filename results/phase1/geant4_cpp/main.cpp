#include <G4Box.hh>
#include <G4Element.hh>
#include <G4EmCalculator.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4Event.hh>
#include <G4Gamma.hh>
#include <G4LogicalVolume.hh>
#include <G4MTRunManager.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4Neutron.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4PhysicalConstants.hh>
#include <G4PVPlacement.hh>
#include <G4Run.hh>
#include <G4RunManagerFactory.hh>
#include <G4Step.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThermalNeutrons.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4UImanager.hh>
#include <G4VModularPhysicsList.hh>
#include <G4VPhysicalVolume.hh>
#include <G4Version.hh>
#include <G4VUserActionInitialization.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4ios.hh>
#include <Randomize.hh>
#include <Shielding.hh>

#include <algorithm>
#include <atomic>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
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
#include <unordered_set>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
constexpr int kRequiredThreads = 16;
constexpr double kPi = 3.1415926535897932384626433832795;

static std::string trim(std::string s) {
    auto notSpace = [](unsigned char c){ return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

static std::vector<std::string> split(const std::string& s, char delim=',') {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) out.push_back(trim(item));
    return out;
}

static std::map<std::string,std::string> readKv(const fs::path& p) {
    std::ifstream in(p);
    if (!in) throw std::runtime_error("Cannot open config: " + p.string());
    std::map<std::string,std::string> kv;
    std::string line;
    while (std::getline(in,line)) {
        line=trim(line);
        if (line.empty() || line[0]=='#') continue;
        auto pos=line.find('=');
        if (pos==std::string::npos) continue;
        kv[trim(line.substr(0,pos))]=trim(line.substr(pos+1));
    }
    return kv;
}

static std::string req(const std::map<std::string,std::string>& kv,const std::string& k) {
    auto it=kv.find(k); if(it==kv.end()) throw std::runtime_error("Missing config key: "+k); return it->second;
}
static double getd(const std::map<std::string,std::string>& kv,const std::string& k,double d=0) {
    auto it=kv.find(k); return it==kv.end()?d:std::stod(it->second);
}
static long long geti64(const std::map<std::string,std::string>& kv,const std::string& k,long long d=0) {
    auto it=kv.find(k); return it==kv.end()?d:std::stoll(it->second);
}

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    std::map<std::string,size_t> idx;
};
static CsvTable readCsv(const fs::path& p) {
    std::ifstream in(p);
    if(!in) throw std::runtime_error("Cannot open CSV: "+p.string());
    CsvTable t; std::string line;
    if(!std::getline(in,line)) throw std::runtime_error("Empty CSV: "+p.string());
    t.header=split(line);
    for(size_t i=0;i<t.header.size();++i) t.idx[t.header[i]]=i;
    while(std::getline(in,line)) {
        if(trim(line).empty()) continue;
        auto r=split(line);
        if(r.size()<t.header.size()) r.resize(t.header.size());
        t.rows.push_back(std::move(r));
    }
    return t;
}
static std::string csvGet(const CsvTable& t,const std::vector<std::string>& r,const std::string& c) {
    auto it=t.idx.find(c); if(it==t.idx.end()) throw std::runtime_error("CSV missing column: "+c); return r.at(it->second);
}

struct SourceBin { double low=0, high=0, prob=0, integral=0; };
struct SourceSpectrum {
    std::vector<SourceBin> bins;
    std::vector<double> cdf;
    double totalIntegral=0;
    static SourceSpectrum load(const fs::path& p) {
        CsvTable t=readCsv(p); SourceSpectrum s; double cum=0;
        for(const auto& r:t.rows) {
            SourceBin b;
            b.low=std::stod(csvGet(t,r,"energy_low_MeV"));
            b.high=std::stod(csvGet(t,r,"energy_high_MeV"));
            b.prob=std::stod(csvGet(t,r,"sampling_probability"));
            b.integral=std::stod(csvGet(t,r,"relative_bin_integral"));
            if(!(b.high>b.low) || b.prob<0) throw std::runtime_error("Invalid source bin");
            s.totalIntegral += b.integral;
            cum += b.prob; s.bins.push_back(b); s.cdf.push_back(cum);
        }
        if(s.bins.empty() || !(cum>0)) throw std::runtime_error("Invalid source spectrum");
        for(auto& x:s.cdf) x/=cum;
        return s;
    }
    double sample() const {
        double u=G4UniformRand(); auto it=std::lower_bound(cdf.begin(),cdf.end(),u);
        size_t i=std::min<size_t>(std::distance(cdf.begin(),it),bins.size()-1);
        const auto& b=bins[i];
        double e=b.low+(b.high-b.low)*G4UniformRand();
        return std::max(e,1.0e-11);
    }
};

struct EnergyBin { double low=0, high=0, du=0; };
static std::vector<EnergyBin> loadBins(const fs::path& p) {
    CsvTable t=readCsv(p); std::vector<EnergyBin> out;
    for(const auto& r:t.rows) {
        EnergyBin b; b.low=std::stod(csvGet(t,r,"energy_lower_MeV")); b.high=std::stod(csvGet(t,r,"energy_upper_MeV"));
        if(!(b.high>b.low) || !(b.low>0)) throw std::runtime_error("Energy-bin bounds must satisfy 0<low<high");
        b.du=std::log(b.high/b.low); out.push_back(b);
    }
    return out;
}

struct Config {
    std::string runId, resultType;
    int sourceProtonMeV=0;
    double shieldCm=0, extraIronCm=0, sourceDistanceBaseCm=400, apertureRadiusCm=5.45;
    double detectorRadiusCm=6.35, detectorLengthCm=12.7, detectorGapCm=0.001;
    std::string transmissionScoringModel;
    std::vector<double> offAxisCm;
    double sourceIntensity=0, peakLow=0, peakHigh=0;
    long long histories=0, seed=1;
    fs::path sourceCsv,binsCsv,outputCsv,normOutputCsv;
    static Config load(const fs::path& p) {
        auto kv=readKv(p); Config c;
        c.runId=req(kv,"run_id"); c.resultType=req(kv,"result_type");
        c.sourceProtonMeV=(int)geti64(kv,"source_proton_mev"); c.shieldCm=getd(kv,"shield_thickness_cm");
        c.extraIronCm=getd(kv,"extra_iron_cm"); c.sourceDistanceBaseCm=getd(kv,"source_distance_base_cm",400);
        c.apertureRadiusCm=getd(kv,"aperture_radius_cm",5.45); c.detectorRadiusCm=getd(kv,"detector_radius_cm",6.35);
        c.detectorLengthCm=getd(kv,"detector_length_cm",12.7); c.detectorGapCm=getd(kv,"detector_gap_cm",0.001);
        if(c.resultType=="neutron_transmission") {
            c.transmissionScoringModel=req(kv,"transmission_scoring_model");
            if(c.transmissionScoringModel!="SINBAD_BC501A_CYLINDRICAL_TRACK_LENGTH_V1")
                throw std::runtime_error("Unsupported neutron-transmission scoring model: "+c.transmissionScoringModel);
            if(!(c.detectorRadiusCm>0.0 && c.detectorLengthCm>0.0 && c.detectorGapCm>=0.0))
                throw std::runtime_error("Invalid cylindrical flux-tally geometry");
        }
        c.sourceIntensity=getd(kv,"source_intensity_n_sr_uC");
        c.peakLow=getd(kv,"peak_low_mev"); c.peakHigh=getd(kv,"peak_high_mev"); c.histories=geti64(kv,"histories"); c.seed=geti64(kv,"random_seed",1);
        for(const auto& x:split(req(kv,"off_axis_cm"))) if(!x.empty()) c.offAxisCm.push_back(std::stod(x));
        if(c.offAxisCm.empty()) c.offAxisCm={0.0};
        c.sourceCsv=req(kv,"source_csv");
        auto it=kv.find("bins_csv"); if(it!=kv.end() && !it->second.empty() && it->second!="NONE") c.binsCsv=it->second;
        c.outputCsv=req(kv,"output_csv"); c.normOutputCsv=req(kv,"normalization_output_csv");
        if(c.histories<1000000) throw std::runtime_error("Phase I stochastic run requires >=1,000,000 histories");
        return c;
    }
};

static G4Material* makeJaeriConcrete() {
    auto* n=G4NistManager::Instance();
    auto* m=new G4Material("JAERI_TIARA_CONCRETE",2.31*g/cm3,9);
    // Add explicitly; the values below are mass fractions converted from JAERI atomic number densities.
    m->AddElement(n->FindOrBuildElement("H"), 0.010859819142947771);
    m->AddElement(n->FindOrBuildElement("O"), 0.48192073243392297);
    m->AddElement(n->FindOrBuildElement("Na"),0.020338354978045932);
    m->AddElement(n->FindOrBuildElement("Mg"),0.010838356046250068);
    m->AddElement(n->FindOrBuildElement("Al"),0.060547665443088476);
    m->AddElement(n->FindOrBuildElement("Si"),0.2242235568799872);
    m->AddElement(n->FindOrBuildElement("K"), 0.010686059058416074);
    m->AddElement(n->FindOrBuildElement("Ca"),0.12395115996130653);
    m->AddElement(n->FindOrBuildElement("Fe"),0.056634296056035024);
    return m;
}

static G4Material* makeBC501A() {
    auto* n=G4NistManager::Instance();
    // BC501A manufacturer data used by the TIARA/SINBAD detector model:
    // density 0.874 g/cm3; H and C atomic densities 0.0482 and 0.0398
    // atom/(barn cm), respectively. Convert those atomic densities to mass
    // fractions while preserving the measured H:C atom ratio.
    constexpr double hNumber=0.0482;
    constexpr double cNumber=0.0398;
    constexpr double hAtomicMass=1.00794;
    constexpr double cAtomicMass=12.0107;
    const double hMass=hNumber*hAtomicMass;
    const double cMass=cNumber*cAtomicMass;
    const double totalMass=hMass+cMass;
    auto* m=new G4Material("BC501A_LIQUID_SCINTILLATOR",0.874*g/cm3,2,kStateLiquid);
    m->AddElement(n->FindOrBuildElement("H"),hMass/totalMass);
    m->AddElement(n->FindOrBuildElement("C"),cMass/totalMass);
    return m;
}

static G4Material* makeNistConcrete() {
    auto* n=G4NistManager::Instance(); auto* m=new G4Material("NIST_ORDINARY_CONCRETE",2.30*g/cm3,10);
    m->AddElement(n->FindOrBuildElement("H"),0.022100); m->AddElement(n->FindOrBuildElement("C"),0.002484);
    m->AddElement(n->FindOrBuildElement("O"),0.574930); m->AddElement(n->FindOrBuildElement("Na"),0.015208);
    m->AddElement(n->FindOrBuildElement("Mg"),0.001266); m->AddElement(n->FindOrBuildElement("Al"),0.019953);
    m->AddElement(n->FindOrBuildElement("Si"),0.304627); m->AddElement(n->FindOrBuildElement("K"),0.010045);
    m->AddElement(n->FindOrBuildElement("Ca"),0.042951); m->AddElement(n->FindOrBuildElement("Fe"),0.006435);
    return m;
}

static std::atomic<long long> gP001GeneratedEvents{0};

class P001Detector final: public G4VUserDetectorConstruction {
    G4Material* material_;
public:
    explicit P001Detector(G4Material* material):material_(material) {
        if(!material_) throw std::runtime_error("P001Detector received null material");
    }
    G4VPhysicalVolume* Construct() override {
        auto* solid=new G4Box("P001_world",1*m,1*m,1*m);
        auto* logical=new G4LogicalVolume(solid,material_,"P001_world");
        return new G4PVPlacement(nullptr,{},logical,"P001_world",nullptr,false,0);
    }
};

class P001PhysicsList final: public G4VModularPhysicsList {
public:
    P001PhysicsList() {
        RegisterPhysics(new G4EmStandardPhysics_option4());
    }
    void SetCuts() override { SetCutsWithDefault(); }
};

class P001PrimaryAction final: public G4VUserPrimaryGeneratorAction {
    G4ParticleGun gun_{1};
public:
    P001PrimaryAction() {
        gun_.SetParticleDefinition(G4Gamma::Gamma());
        gun_.SetParticleEnergy(1.0*MeV);
        gun_.SetParticlePosition(G4ThreeVector(0,0,0));
        gun_.SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    }
    void GeneratePrimaries(G4Event* event) override {
        ++gP001GeneratedEvents;
        gun_.GeneratePrimaryVertex(event);
    }
};

class P001Actions final: public G4VUserActionInitialization {
public:
    void BuildForMaster() const override {}
    void Build() const override {
        SetUserAction(new P001PrimaryAction());
    }
};

class BenchmarkDetector final: public G4VUserDetectorConstruction {
    Config cfg_;
public: explicit BenchmarkDetector(Config c):cfg_(std::move(c)){}
    G4VPhysicalVolume* Construct() override {
        auto* n=G4NistManager::Instance(); auto* vac=n->FindOrBuildMaterial("G4_Galactic"); auto* iron=n->FindOrBuildMaterial("G4_Fe"); auto* concrete=makeJaeriConcrete(); auto* bc501a=makeBC501A();
        const double sourceZcm=-(cfg_.sourceDistanceBaseCm+cfg_.extraIronCm);
        const double zMax=cfg_.shieldCm+cfg_.detectorGapCm+100.0;
        const double halfZ=std::max(std::abs(sourceZcm),zMax)+100.0;
        auto* ws=new G4Box("world",400*cm,400*cm,halfZ*cm); auto* wl=new G4LogicalVolume(ws,vac,"world"); auto* wp=new G4PVPlacement(nullptr,{},wl,"world",nullptr,false,0);
        if(cfg_.extraIronCm>0) {
            auto* outer=new G4Box("iron_outer",60*cm,60*cm,0.5*cfg_.extraIronCm*cm);
            auto* hole=new G4Tubs("iron_hole",0,cfg_.apertureRadiusCm*cm,0.6*cfg_.extraIronCm*cm,0,twopi);
            auto* solid=new G4SubtractionSolid("iron_collimator",outer,hole);
            auto* logical=new G4LogicalVolume(solid,iron,"iron_collimator");
            new G4PVPlacement(nullptr,G4ThreeVector(0,0,-0.5*cfg_.extraIronCm*cm),logical,"iron_collimator",wl,false,0);
        }
        if(cfg_.shieldCm>0) {
            auto* cs=new G4Box("concrete",60*cm,60*cm,0.5*cfg_.shieldCm*cm); auto* cl=new G4LogicalVolume(cs,concrete,"concrete");
            new G4PVPlacement(nullptr,G4ThreeVector(0,0,0.5*cfg_.shieldCm*cm),cl,"concrete",wl,false,0);
        }
        if(cfg_.resultType=="neutron_transmission") {
            // SINBAD/TIARA calculation model: cylindrical flux estimators matching
            // the 12.7-cm-diameter x 12.7-cm-long BC501A measurement geometry.
            // G4Tubs is aligned with z by construction.  The detector is filled
            // with explicit BC501A liquid scintillator, matching the SINBAD F4
            // volume-tally model; neutron track length is accumulated in this volume.
            auto* ds=new G4Tubs("jaeri_bc501a_flux_tally_solid",0,cfg_.detectorRadiusCm*cm,
                                0.5*cfg_.detectorLengthCm*cm,0,twopi);
            auto* dl=new G4LogicalVolume(ds,bc501a,"jaeri_bc501a_flux_tally_logical");
            const double detectorCenterZcm=cfg_.shieldCm+cfg_.detectorGapCm+0.5*cfg_.detectorLengthCm;
            for(size_t d=0; d<cfg_.offAxisCm.size(); ++d) {
                new G4PVPlacement(nullptr,
                    G4ThreeVector(cfg_.offAxisCm[d]*cm,0,detectorCenterZcm*cm),
                    dl,"jaeri_bc501a_flux_tally",wl,false,static_cast<int>(d));
            }
        }
        return wp;
    }
};

struct Accum { double sum=0,sumsq=0; long long n=0; void add(double x){sum+=x;sumsq+=x*x;++n;} };
struct SharedTallies {
    std::mutex mu; std::vector<Accum> bins; Accum norm,dose;
    size_t nDet=0,nBins=0;
    SharedTallies(size_t nd,size_t nb):bins(nd*nb),nDet(nd),nBins(nb){}
    void merge(const std::vector<Accum>& b,const Accum& no,const Accum& do_) { std::lock_guard<std::mutex> lock(mu); for(size_t i=0;i<bins.size();++i){bins[i].sum+=b[i].sum;bins[i].sumsq+=b[i].sumsq;bins[i].n+=b[i].n;} norm.sum+=no.sum;norm.sumsq+=no.sumsq;norm.n+=no.n; dose.sum+=do_.sum;dose.sumsq+=do_.sumsq;dose.n+=do_.n; }
};

class WorkerRunAction;
class WorkerEventAction final: public G4UserEventAction {
    WorkerRunAction* run_; std::vector<double> eventBins_; std::vector<size_t> touchedBins_; double eventNorm_=0,eventDose_=0;
    std::unordered_set<int> outputTracks_,normTracks_;
public:
    WorkerEventAction(WorkerRunAction* r,size_t n):run_(r),eventBins_(n,0){touchedBins_.reserve(std::min<size_t>(n,16));}
    void BeginOfEventAction(const G4Event*) override {
        // Exact sparse reset: bins not touched by the previous event are already zero.
        for(size_t i:touchedBins_) eventBins_[i]=0.0;
        touchedBins_.clear(); eventNorm_=0; eventDose_=0; outputTracks_.clear(); normTracks_.clear();
    }
    void EndOfEventAction(const G4Event*) override;
    void addBin(size_t i,double x){
        if(!(x>0.0)) return;
        double& v=eventBins_.at(i);
        if(v==0.0) touchedBins_.push_back(i);
        v+=x;
    }
    void addNorm(double x){eventNorm_+=x;} void addDose(double x){eventDose_+=x;}
    bool firstOutput(int id){return outputTracks_.insert(id).second;} bool firstNorm(int id){return normTracks_.insert(id).second;}
    const std::vector<size_t>& touchedBins() const {return touchedBins_;}
};

class WorkerRunAction final: public G4UserRunAction {
    SharedTallies& shared_; std::vector<Accum> bins_; Accum norm_,dose_; long long eventCount_=0;
public: WorkerRunAction(SharedTallies& s):shared_(s),bins_(s.bins.size()){}
    void addEvent(const std::vector<double>& b,const std::vector<size_t>& touched,double no,double do_) {
        // Untouched bins are exact zero scores. Omitting zero additions leaves sum and sumsq
        // bit-for-bit unchanged; the full event count is restored once at end-of-run.
        for(size_t i:touched) bins_.at(i).add(b.at(i));
        ++eventCount_; norm_.add(no); dose_.add(do_);
    }
    void EndOfRunAction(const G4Run*) override {
        for(auto& a:bins_) a.n=eventCount_;
        shared_.merge(bins_,norm_,dose_);
    }
};
void WorkerEventAction::EndOfEventAction(const G4Event*) {run_->addEvent(eventBins_,touchedBins_,eventNorm_,eventDose_);} 

static double icrp21Coeff(double e) {
    static const std::array<double,16> E={2.5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,0.5,1,2,5,10,20,50,100};
    static const std::array<double,16> H={10.68,11.57,12.63,12.08,11.57,10.29,9.92,57.87,198.41,326.80,396.83,408.50,408.50,427.35,455.37,496.03};
    if(e<=E.front()) return H.front(); if(e>=E.back()) return H.back();
    auto it=std::upper_bound(E.begin(),E.end(),e); size_t j=std::distance(E.begin(),it), i=j-1;
    double x=(std::log(e)-std::log(E[i]))/(std::log(E[j])-std::log(E[i]));
    return std::exp(std::log(H[i])+x*(std::log(H[j])-std::log(H[i])));
}

class PrimaryAction final: public G4VUserPrimaryGeneratorAction {
    Config cfg_; SourceSpectrum src_; G4ParticleGun gun_{1}; double distanceCm_;
public: PrimaryAction(Config c,SourceSpectrum s):cfg_(std::move(c)),src_(std::move(s)),distanceCm_(cfg_.sourceDistanceBaseCm+cfg_.extraIronCm) {
    gun_.SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("neutron"));
    gun_.SetParticlePosition(G4ThreeVector(0,0,-distanceCm_*cm));
    }
    void GeneratePrimaries(G4Event* evt) override {
        double e=src_.sample(); gun_.SetParticleEnergy(e*MeV);
        const double thetaMax=std::atan(cfg_.apertureRadiusCm/cfg_.sourceDistanceBaseCm); const double cosMax=std::cos(thetaMax);
        double cosT=1.0-G4UniformRand()*(1.0-cosMax); double sinT=std::sqrt(std::max(0.0,1.0-cosT*cosT)); double phi=twopi*G4UniformRand();
        gun_.SetParticleMomentumDirection(G4ThreeVector(sinT*std::cos(phi),sinT*std::sin(phi),cosT)); gun_.GeneratePrimaryVertex(evt);
    }
};

class StepAction final: public G4UserSteppingAction {
    Config cfg_; const std::vector<EnergyBin>& bins_; WorkerEventAction& event_;
    double normZ_=0.0,outZ_=0.0,apertureRadiusSq_=0.0,detectorRadiusSqCm2_=0.0,detectorAreaCm2_=0.0,detectorVolumeCm3_=0.0,normAreaCm2_=0.0;
    bool binsOrderedNonoverlapping_=true;
public:
    StepAction(Config c,const std::vector<EnergyBin>& b,WorkerEventAction& e):cfg_(std::move(c)),bins_(b),event_(e) {
        normZ_=-0.001*cm;
        outZ_=(cfg_.shieldCm+cfg_.detectorGapCm)*cm;
        apertureRadiusSq_=std::pow(cfg_.apertureRadiusCm*cm,2);
        detectorRadiusSqCm2_=cfg_.detectorRadiusCm*cfg_.detectorRadiusCm;
        detectorAreaCm2_=kPi*detectorRadiusSqCm2_;
        detectorVolumeCm3_=detectorAreaCm2_*cfg_.detectorLengthCm;
        normAreaCm2_=kPi*cfg_.apertureRadiusCm*cfg_.apertureRadiusCm;
        for(size_t i=1;i<bins_.size();++i) {
            if(bins_[i].low<bins_[i-1].high || bins_[i].high<bins_[i-1].high) {
                binsOrderedNonoverlapping_=false; break;
            }
        }
    }
    static bool crosses(const G4Step* st,double z,G4ThreeVector& pos,double& energy,double& cosz) {
        auto* pre=st->GetPreStepPoint(); auto* post=st->GetPostStepPoint(); double z0=pre->GetPosition().z(), z1=post->GetPosition().z();
        if(!(z0<z && z1>=z)) return false; double f=(z-z0)/(z1-z0); pos=pre->GetPosition()+f*(post->GetPosition()-pre->GetPosition());
        auto dir=pre->GetMomentumDirection(); cosz=dir.z(); if(!(cosz>1e-9)) return false; energy=pre->GetKineticEnergy()/MeV; return true;
    }
    int energyBin(double e) const {
        // Fast path for the canonical ordered non-overlapping benchmark bins.
        // The fallback preserves the former first-match linear semantics for any
        // unexpected legacy bin ordering or overlap.
        if(!binsOrderedNonoverlapping_) {
            for(size_t k=0;k<bins_.size();++k) if(e>=bins_[k].low && e<bins_[k].high) return static_cast<int>(k);
            return -1;
        }
        size_t lo=0,hi=bins_.size();
        while(lo<hi){size_t mid=lo+(hi-lo)/2; if(e<bins_[mid].high) hi=mid; else lo=mid+1;}
        if(lo<bins_.size() && e>=bins_[lo].low && e<bins_[lo].high) return static_cast<int>(lo);
        return -1;
    }
    void UserSteppingAction(const G4Step* st) override {
        if(st->GetTrack()->GetDefinition()!=G4Neutron::NeutronDefinition()) return;
        G4ThreeVector pos; double e=0,cosz=0; int tid=st->GetTrack()->GetTrackID();
        if(crosses(st,normZ_,pos,e,cosz) && event_.firstNorm(tid)) {
            double r2=pos.x()*pos.x()+pos.y()*pos.y(); if(r2<=apertureRadiusSq_ && e>=cfg_.peakLow && e<=cfg_.peakHigh) {
                event_.addNorm(1.0/(normAreaCm2_*cosz));
            }
        }
        if(cfg_.resultType=="neutron_transmission") {
            // SINBAD/TIARA scalar-fluence estimator: sum neutron path length in
            // the cylindrical BC501A tally volume and divide by its volume.
            // No direction cut and no first-crossing suppression are applied.
            auto* pre=st->GetPreStepPoint();
            auto* pv=pre->GetTouchableHandle()->GetVolume();
            if(!pv || pv->GetName()!="jaeri_bc501a_flux_tally") return;
            const int d=pv->GetCopyNo();
            if(d<0 || static_cast<size_t>(d)>=cfg_.offAxisCm.size()) return;
            const double stepLengthCm=st->GetStepLength()/cm;
            if(!(stepLengthCm>0.0)) return;
            e=pre->GetKineticEnergy()/MeV;
            const int k=energyBin(e);
            if(k<0) return;
            const double flu=st->GetTrack()->GetWeight()*stepLengthCm/detectorVolumeCm3_;
            event_.addBin(static_cast<size_t>(d)*bins_.size()+static_cast<size_t>(k),
                          flu/bins_[static_cast<size_t>(k)].du);
            return;
        }

        // Legacy plane scorer retained only for non-transmission observables
        // whose v12.9 semantics are intentionally unchanged (ICRP-21 and the
        // controlled fixed-geometry discovery sweep).
        if(!crosses(st,outZ_,pos,e,cosz) || !event_.firstOutput(tid)) return;
        const int k=cfg_.resultType=="icrp21_dose" ? -1 : energyBin(e);
        for(size_t d=0;d<cfg_.offAxisCm.size();++d) {
            double dx=pos.x()/cm-cfg_.offAxisCm[d], dy=pos.y()/cm; if(dx*dx+dy*dy>detectorRadiusSqCm2_) continue;
            double flu=1.0/(detectorAreaCm2_*cosz);
            if(cfg_.resultType=="icrp21_dose") { if(d==0) event_.addDose(flu*icrp21Coeff(e)); }
            else if(k>=0) event_.addBin(d*bins_.size()+static_cast<size_t>(k),flu/bins_[static_cast<size_t>(k)].du);
        }
    }
};

class Actions final: public G4VUserActionInitialization {
    Config cfg_; SourceSpectrum src_; std::vector<EnergyBin> bins_; SharedTallies& shared_;
public: Actions(Config c,SourceSpectrum s,std::vector<EnergyBin> b,SharedTallies& sh):cfg_(std::move(c)),src_(std::move(s)),bins_(std::move(b)),shared_(sh){}
    void BuildForMaster() const override {SetUserAction(new WorkerRunAction(shared_));}
    void Build() const override {auto* r=new WorkerRunAction(shared_); auto* e=new WorkerEventAction(r,shared_.bins.size()); SetUserAction(new PrimaryAction(cfg_,src_));SetUserAction(r);SetUserAction(e);SetUserAction(new StepAction(cfg_,bins_,*e));}
};

static std::pair<double,double> meanSem(const Accum& a,long long N) {
    if(N<=0) return {0,0}; double mean=a.sum/N; if(N<=1) return {mean,0}; double var=(a.sumsq-N*mean*mean)/(N-1); var=std::max(0.0,var); return {mean,std::sqrt(var/N)};
}

static void writeNeutronOutputs(const Config& cfg,const SourceSpectrum& src,const std::vector<EnergyBin>& bins,const SharedTallies& sh) {
    double sourceConeDistance=cfg.sourceDistanceBaseCm; double theta=std::atan(cfg.apertureRadiusCm/sourceConeDistance); double omega=2*kPi*(1-std::cos(theta));
    double physicalSourceNeutronsPerUC=cfg.sourceIntensity*src.totalIntegral*omega; double primaryWeightPerUC=physicalSourceNeutronsPerUC/cfg.histories;
    fs::create_directories(cfg.outputCsv.parent_path());
    std::ofstream out(cfg.outputCsv); out<<std::setprecision(17);
    if(cfg.resultType=="icrp21_dose") {
        auto [m,s]=meanSem(sh.dose,cfg.histories); out<<"run_id,shield_thickness_cm,mc_icrp21_dose_equivalent_uSv_per_uC,mc_sigma_icrp21_dose_equivalent_uSv_per_uC\n";
        out<<cfg.runId<<","<<cfg.shieldCm<<","<<m*physicalSourceNeutronsPerUC*1e-6<<","<<s*physicalSourceNeutronsPerUC*1e-6<<"\n";
    } else {
        out<<"run_id,off_axis_cm,energy_lower_MeV,energy_upper_MeV,mc_lethargy_flux_n_cm2_per_uC,mc_sigma_lethargy_flux_n_cm2_per_uC\n";
        for(size_t d=0;d<cfg.offAxisCm.size();++d) for(size_t k=0;k<bins.size();++k) {auto [m,s]=meanSem(sh.bins[d*bins.size()+k],cfg.histories); out<<cfg.runId<<","<<cfg.offAxisCm[d]<<","<<bins[k].low<<","<<bins[k].high<<","<<m*physicalSourceNeutronsPerUC<<","<<s*physicalSourceNeutronsPerUC<<"\n";}
    }
    std::ofstream no(cfg.normOutputCsv); no<<std::setprecision(17); auto [nm,ns]=meanSem(sh.norm,cfg.histories);
    no<<"run_id,mc_incident_peak_fluence_n_cm2_per_uC,mc_sigma_incident_peak_fluence_n_cm2_per_uC,source_solid_angle_sr,normalization_scale_per_primary_per_uC,normalization_multiplier_for_per_primary_mean_per_uC\n";
    no<<cfg.runId<<","<<nm*physicalSourceNeutronsPerUC<<","<<ns*physicalSourceNeutronsPerUC<<","<<omega<<","<<primaryWeightPerUC<<","<<physicalSourceNeutronsPerUC<<"\n";
}


struct DepthBin { double center=0, low=0, high=0; };
static std::vector<DepthBin> loadDepthBinsFromReference(const fs::path& p) {
    CsvTable t=readCsv(p); std::vector<double> centers;
    for(const auto& r:t.rows) centers.push_back(std::stod(csvGet(t,r,"depth_cm")));
    std::sort(centers.begin(),centers.end());
    centers.erase(std::unique(centers.begin(),centers.end()),centers.end());
    if(centers.size()<2 || centers.front()<0) throw std::runtime_error("PDD reference requires >=2 nonnegative depth coordinates");
    std::vector<double> edges(centers.size()+1,0.0);
    edges[0]=std::max(0.0,centers[0]-0.5*(centers[1]-centers[0]));
    for(size_t i=1;i<centers.size();++i) edges[i]=0.5*(centers[i-1]+centers[i]);
    edges.back()=centers.back()+0.5*(centers.back()-centers[centers.size()-2]);
    std::vector<DepthBin> out; out.reserve(centers.size());
    for(size_t i=0;i<centers.size();++i) {
        if(!(edges[i+1]>edges[i])) throw std::runtime_error("Invalid PDD depth bin width");
        out.push_back({centers[i],edges[i],edges[i+1]});
    }
    return out;
}


struct PddConfig {
    std::string runId,sourceModelId,geometryMatchStatus,sourcePlaneMode,phaseSpaceModel;
    double ssdCm=100,fieldXCm=40,fieldYCm=40,scoreHalfWidthCm=1,waterHalfXYCm=40,waterDepthCm=50,normalizationDepthCm=0,sourcePlaneOffsetCm=0.001;
    bool phaseSpaceComplete=false;
    long long histories=0,seed=1;
    fs::path sourceCsv,referenceCsv,outputCsv;
    static PddConfig load(const fs::path& p) {
        auto kv=readKv(p); PddConfig c;
        c.runId=req(kv,"run_id"); c.sourceModelId=req(kv,"source_model_id");
        c.ssdCm=getd(kv,"ssd_cm",100); c.fieldXCm=getd(kv,"field_x_cm",40); c.fieldYCm=getd(kv,"field_y_cm",40);
        c.scoreHalfWidthCm=getd(kv,"score_half_width_cm",1); c.waterHalfXYCm=getd(kv,"water_half_xy_cm",40); c.waterDepthCm=getd(kv,"water_depth_cm",50);
        c.normalizationDepthCm=getd(kv,"normalization_depth_cm"); c.sourcePlaneOffsetCm=getd(kv,"source_plane_offset_cm",0.001);
        c.histories=geti64(kv,"histories"); c.seed=geti64(kv,"random_seed",1);
        c.sourceCsv=req(kv,"source_csv"); c.referenceCsv=req(kv,"reference_csv"); c.outputCsv=req(kv,"output_csv");
        c.sourcePlaneMode=req(kv,"source_plane_mode"); c.phaseSpaceModel=req(kv,"phase_space_model");
        c.phaseSpaceComplete=(geti64(kv,"phase_space_complete",0)!=0);
        auto it=kv.find("geometry_match_status"); if(it!=kv.end()) c.geometryMatchStatus=it->second;
        if(c.histories<1000000) throw std::runtime_error("Phase I photon PDD run requires >=1,000,000 histories");
        if(!(c.ssdCm>0 && c.fieldXCm>0 && c.fieldYCm>0 && c.scoreHalfWidthCm>0 && c.waterDepthCm>0 && c.sourcePlaneOffsetCm>0)) throw std::runtime_error("Invalid photon PDD geometry");
        if(c.sourcePlaneMode!="water_surface_factorized") throw std::runtime_error("Photon PDD source must be launched from the water-surface phase-space plane");
        if(c.phaseSpaceModel!="aggregate_energy_uniform_xy_virtual_source_divergence") throw std::runtime_error("Unsupported photon PDD phase-space factorization");
        if(c.phaseSpaceComplete) throw std::runtime_error("Current 1-D production spectrum does not contain complete x/y/energy/angle phase-space correlations");
        return c;
    }
};

class PddDetector final: public G4VUserDetectorConstruction {
    PddConfig cfg_;
public: explicit PddDetector(PddConfig c):cfg_(std::move(c)){}
    G4VPhysicalVolume* Construct() override {
        auto* n=G4NistManager::Instance(); auto* air=n->FindOrBuildMaterial("G4_AIR"); auto* water=n->FindOrBuildMaterial("G4_WATER");
        // Surface-phase-space contract:
        //   water surface z=0;
        //   primary vertices are placed a small numerical offset in air at z=-sourcePlaneOffset;
        //   SSD is used only as a virtual-source distance to assign divergence.
        // No extra 100-cm air transport is applied to a spectrum already defined at the phantom surface.
        const double worldHalfXY=std::max(100.0,cfg_.waterHalfXYCm+10.0);
        const double worldHalfZ=std::max(25.0,cfg_.waterDepthCm+20.0);
        auto* ws=new G4Box("pdd_world",worldHalfXY*cm,worldHalfXY*cm,worldHalfZ*cm); auto* wl=new G4LogicalVolume(ws,air,"pdd_world");
        auto* wp=new G4PVPlacement(nullptr,G4ThreeVector(),wl,"pdd_world",nullptr,false,0);
        auto* ps=new G4Box("pdd_water",cfg_.waterHalfXYCm*cm,cfg_.waterHalfXYCm*cm,0.5*cfg_.waterDepthCm*cm);
        auto* pl=new G4LogicalVolume(ps,water,"pdd_water");
        new G4PVPlacement(nullptr,G4ThreeVector(0,0,0.5*cfg_.waterDepthCm*cm),pl,"pdd_water",wl,false,0);
        return wp;
    }
};

struct PddSharedTallies {
    std::mutex mu;
    std::vector<Accum> dose;
    std::vector<double> crossWithNormSum;
    long long crossN=0;
    size_t normIndex=0;
    PddSharedTallies(size_t n,size_t norm):dose(n),crossWithNormSum(n,0.0),normIndex(norm){}
    void merge(const std::vector<Accum>& x,const std::vector<double>& cross,long long n) {
        std::lock_guard<std::mutex> lock(mu);
        for(size_t i=0;i<dose.size();++i){
            dose[i].sum+=x[i].sum;
            dose[i].sumsq+=x[i].sumsq;
            dose[i].n+=x[i].n;
            crossWithNormSum[i]+=cross[i];
        }
        crossN+=n;
    }
};

class PddRunAction;
class PddEventAction final: public G4UserEventAction {
    PddRunAction* run_; std::vector<double> eventDose_;
public:
    PddEventAction(PddRunAction* r,size_t n):run_(r),eventDose_(n,0){}
    void BeginOfEventAction(const G4Event*) override {std::fill(eventDose_.begin(),eventDose_.end(),0.0);}
    void EndOfEventAction(const G4Event*) override;
    void addDose(size_t i,double x){eventDose_.at(i)+=x;}
};

class PddRunAction final: public G4UserRunAction {
    PddSharedTallies& shared_;
    std::vector<Accum> dose_;
    std::vector<double> crossWithNormSum_;
    long long crossN_=0;
public:
    explicit PddRunAction(PddSharedTallies& s):shared_(s),dose_(s.dose.size()),crossWithNormSum_(s.dose.size(),0.0){}
    void addEvent(const std::vector<double>& x){
        const double yn=x.at(shared_.normIndex);
        for(size_t i=0;i<x.size();++i){
            dose_[i].add(x[i]);
            crossWithNormSum_[i]+=x[i]*yn;
        }
        ++crossN_;
    }
    void EndOfRunAction(const G4Run*) override {shared_.merge(dose_,crossWithNormSum_,crossN_);}
};
void PddEventAction::EndOfEventAction(const G4Event*) {run_->addEvent(eventDose_);}

class PddPrimaryAction final: public G4VUserPrimaryGeneratorAction {
    PddConfig cfg_; SourceSpectrum src_; G4ParticleGun gun_{1};
public:
    PddPrimaryAction(PddConfig c,SourceSpectrum s):cfg_(std::move(c)),src_(std::move(s)) {
        gun_.SetParticleDefinition(G4Gamma::Gamma());
    }
    void GeneratePrimaries(G4Event* evt) override {
        gun_.SetParticleEnergy(src_.sample()*MeV);
        const double x=(G4UniformRand()-0.5)*cfg_.fieldXCm;
        const double y=(G4UniformRand()-0.5)*cfg_.fieldYCm;

        // The production energy spectrum is a surface-plane energy distribution.
        // Sample a factorized lateral coordinate at the water surface and use the
        // reported SSD only to assign the local virtual-source divergence.
        gun_.SetParticlePosition(G4ThreeVector(x*cm,y*cm,-cfg_.sourcePlaneOffsetCm*cm));
        G4ThreeVector dir(x*cm,y*cm,cfg_.ssdCm*cm);
        dir=dir.unit();
        gun_.SetParticleMomentumDirection(dir);
        gun_.GeneratePrimaryVertex(evt);
    }
};

class PddStepAction final: public G4UserSteppingAction {
    PddConfig cfg_; const std::vector<DepthBin>& bins_; PddEventAction& event_;
public:
    PddStepAction(PddConfig c,const std::vector<DepthBin>& b,PddEventAction& e):cfg_(std::move(c)),bins_(b),event_(e){}
    void UserSteppingAction(const G4Step* st) override {
        const double edep=st->GetTotalEnergyDeposit()/MeV; if(!(edep>0)) return;
        auto* pre=st->GetPreStepPoint(); auto touch=pre->GetTouchableHandle(); auto* volume=touch->GetVolume(); if(!volume) return;
        if(volume->GetLogicalVolume()->GetName()!="pdd_water") return;
        auto pos=0.5*(pre->GetPosition()+st->GetPostStepPoint()->GetPosition());
        const double x=pos.x()/cm,y=pos.y()/cm,z=pos.z()/cm;
        if(std::abs(x)>cfg_.scoreHalfWidthCm || std::abs(y)>cfg_.scoreHalfWidthCm || z<0) return;
        for(size_t i=0;i<bins_.size();++i) {
            const bool inside=(z>=bins_[i].low && (z<bins_[i].high || (i+1==bins_.size() && z<=bins_[i].high)));
            if(!inside) continue;
            const double volumeCm3=4.0*cfg_.scoreHalfWidthCm*cfg_.scoreHalfWidthCm*(bins_[i].high-bins_[i].low);
            if(volumeCm3>0) event_.addDose(i,edep/volumeCm3); // relative water dose proxy; density cancels in PDD ratio
            break;
        }
    }
};

class PddActions final: public G4VUserActionInitialization {
    PddConfig cfg_; SourceSpectrum src_; std::vector<DepthBin> bins_; PddSharedTallies& shared_;
public:
    PddActions(PddConfig c,SourceSpectrum s,std::vector<DepthBin> b,PddSharedTallies& sh):cfg_(std::move(c)),src_(std::move(s)),bins_(std::move(b)),shared_(sh){}
    void BuildForMaster() const override {SetUserAction(new PddRunAction(shared_));}
    void Build() const override {
        auto* r=new PddRunAction(shared_);
        auto* e=new PddEventAction(r,shared_.dose.size());
        SetUserAction(new PddPrimaryAction(cfg_,src_));
        SetUserAction(r);
        SetUserAction(e);
        SetUserAction(new PddStepAction(cfg_,bins_,*e));
    }
};

static int runPhotonPdd(const fs::path& configPath) {
#ifndef G4MULTITHREADED
    throw std::runtime_error("This Phase I executable requires a multithreaded Geant4 build; G4MULTITHREADED is not defined.");
#else
    PddConfig cfg=PddConfig::load(configPath);
    SourceSpectrum src=SourceSpectrum::load(cfg.sourceCsv);
    auto bins=loadDepthBinsFromReference(cfg.referenceCsv);

    size_t norm=0;
    double best=std::numeric_limits<double>::infinity();
    for(size_t i=0;i<bins.size();++i){
        double d=std::abs(bins[i].center-cfg.normalizationDepthCm);
        if(d<best){best=d;norm=i;}
    }

    PddSharedTallies shared(bins.size(),norm);
    G4Random::setTheSeed((long)cfg.seed);
    auto* rm=G4RunManagerFactory::CreateRunManager(G4RunManagerType::MTOnly);
    auto* mt=dynamic_cast<G4MTRunManager*>(rm);
    if(!mt) throw std::runtime_error("Geant4 MT run manager unavailable");
    mt->SetNumberOfThreads(kRequiredThreads);
    rm->SetUserInitialization(new PddDetector(cfg));
    class LocalPddPhysics final: public G4VModularPhysicsList {
    public:
        LocalPddPhysics(){SetVerboseLevel(1);RegisterPhysics(new G4EmStandardPhysics_option4());}
    };
    rm->SetUserInitialization(new LocalPddPhysics());
    rm->SetUserInitialization(new PddActions(cfg,src,bins,shared));
    auto* ui=G4UImanager::GetUIpointer();
    ui->ApplyCommand("/control/verbose 1");
    ui->ApplyCommand("/run/verbose 1");
    ui->ApplyCommand("/event/verbose 0");
    ui->ApplyCommand("/tracking/verbose 0");
    G4cout<<"Phase I native C++ photon PDD surface-plane diagnostic runner; worker threads="<<kRequiredThreads<<G4endl;
    rm->Initialize();
    rm->BeamOn((G4int)cfg.histories);

    if(shared.crossN!=cfg.histories) {
        delete rm;
        throw std::runtime_error("Photon PDD event-level covariance tally count does not equal requested histories");
    }

    std::vector<double> mean(bins.size(),0),sem(bins.size(),0);
    for(size_t i=0;i<bins.size();++i){
        auto ms=meanSem(shared.dose[i],cfg.histories);
        mean[i]=ms.first;
        sem[i]=ms.second;
    }

    if(!(mean[norm]>0)) {
        delete rm;
        throw std::runtime_error("Photon PDD normalization bin has zero dose");
    }

    fs::create_directories(cfg.outputCsv.parent_path());
    std::ofstream out(cfg.outputCsv);
    out<<std::setprecision(17);
    out<<"run_id,depth_cm,mc_pdd_percent,mc_sigma_pdd_percent,mc_raw_mean_relative_dose,mc_raw_sem_relative_dose,mc_covariance_with_normalization_mean,normalization_depth_cm,score_half_width_cm,virtual_source_distance_cm,field_size_x_cm,field_size_y_cm,source_plane_mode,source_plane_offset_cm,phase_space_model,phase_space_complete,mc_uncertainty_method\n";

    const double N=static_cast<double>(cfg.histories);
    for(size_t i=0;i<bins.size();++i){
        const double p=100.0*mean[i]/mean[norm];
        double covMean=0.0;
        if(cfg.histories>1){
            const double sampleCov=(shared.crossWithNormSum[i]-N*mean[i]*mean[norm])/(N-1.0);
            covMean=sampleCov/N;
        }

        double sp=0.0;
        if(i!=norm){
            const double dA=100.0/mean[norm];
            const double dB=-100.0*mean[i]/(mean[norm]*mean[norm]);
            const double varA=sem[i]*sem[i];
            const double varB=sem[norm]*sem[norm];
            const double varP=dA*dA*varA+dB*dB*varB+2.0*dA*dB*covMean;
            sp=std::sqrt(std::max(0.0,varP));
        }

        out<<cfg.runId<<","<<bins[i].center<<","<<p<<","<<sp<<","<<mean[i]<<","<<sem[i]<<","<<covMean<<","<<bins[norm].center<<","<<cfg.scoreHalfWidthCm<<","<<cfg.ssdCm<<","<<cfg.fieldXCm<<","<<cfg.fieldYCm<<","<<cfg.sourcePlaneMode<<","<<cfg.sourcePlaneOffsetCm<<","<<cfg.phaseSpaceModel<<","<<(cfg.phaseSpaceComplete?1:0)<<",event_level_delta_method_with_covariance\n";
    }

    out.flush();
    if(!out){
        delete rm;
        throw std::runtime_error("Failed writing photon PDD output");
    }
    delete rm;
    return 0;
#endif
}

static int runNeutron(const fs::path& configPath) {
#ifndef G4MULTITHREADED
    throw std::runtime_error("This Phase I executable requires a multithreaded Geant4 build; G4MULTITHREADED is not defined.");
#else
    Config cfg=Config::load(configPath); SourceSpectrum src=SourceSpectrum::load(cfg.sourceCsv); std::vector<EnergyBin> bins;
    if(cfg.resultType!="icrp21_dose") bins=loadBins(cfg.binsCsv);
    SharedTallies shared(cfg.offAxisCm.size(),bins.size());
    G4Random::setTheSeed((long)cfg.seed);
    auto* rm=G4RunManagerFactory::CreateRunManager(G4RunManagerType::MTOnly); auto* mt=dynamic_cast<G4MTRunManager*>(rm); if(!mt) throw std::runtime_error("Geant4 MT run manager unavailable");
    mt->SetNumberOfThreads(kRequiredThreads);
    rm->SetUserInitialization(new BenchmarkDetector(cfg)); auto* physics=new Shielding(0); physics->RegisterPhysics(new G4ThermalNeutrons()); rm->SetUserInitialization(physics);
    rm->SetUserInitialization(new Actions(cfg,src,bins,shared));
    auto* ui=G4UImanager::GetUIpointer(); ui->ApplyCommand("/control/verbose 1");ui->ApplyCommand("/run/verbose 1");ui->ApplyCommand("/event/verbose 0");ui->ApplyCommand("/tracking/verbose 0");
    G4cout<<"Phase I native C++ neutron runner; worker threads="<<kRequiredThreads<<G4endl;
    rm->Initialize();
    const auto transportStart=std::chrono::steady_clock::now();
    rm->BeamOn((G4int)cfg.histories);
    const double transportSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-transportStart).count();
    const double historiesPerSecond=transportSeconds>0.0 ? static_cast<double>(cfg.histories)/transportSeconds : 0.0;
    G4cout<<"NEUTRON_TRANSPORT_RUNTIME_SECONDS="<<std::setprecision(10)<<transportSeconds<<G4endl;
    G4cout<<"NEUTRON_HISTORIES_PER_SECOND="<<std::setprecision(10)<<historiesPerSecond<<G4endl;
    G4cout<<"NEUTRON_SCORING_IMPLEMENTATION=SPARSE_ZERO_EQUIVALENT_V1"<<G4endl;
    G4cout<<"NEUTRON_SOURCE_CONE_REFERENCE=ROTARY_SHUTTER_EXIT"<<G4endl;
    G4cout<<"NEUTRON_SOURCE_CONE_REFERENCE_DISTANCE_CM="<<std::setprecision(10)<<cfg.sourceDistanceBaseCm<<G4endl;
    G4cout<<"NEUTRON_SOURCE_TO_SHIELD_FRONT_CM="<<std::setprecision(10)<<(cfg.sourceDistanceBaseCm+cfg.extraIronCm)<<G4endl;
    G4cout<<"NEUTRON_ADDITIONAL_IRON_TRANSPORT="<<(cfg.extraIronCm>0.0 ? "EXPLICIT_DOWNSTREAM_HOLLOW_IRON" : "NONE")<<G4endl;
    writeNeutronOutputs(cfg,src,bins,shared); delete rm; return 0;
#endif
}

static int runP001(const fs::path& target,const fs::path& outPath) {
#ifndef G4MULTITHREADED
    throw std::runtime_error("This Phase I executable requires a multithreaded Geant4 build; G4MULTITHREADED is not defined.");
#else
    gP001GeneratedEvents.store(0);

    // Material must exist before EM initialization.
    auto* mat=makeNistConcrete();

    auto* rm=G4RunManagerFactory::CreateRunManager(G4RunManagerType::MTOnly);
    auto* mt=dynamic_cast<G4MTRunManager*>(rm);
    if(!mt) {
        delete rm;
        throw std::runtime_error("Geant4 MT run manager unavailable for P001");
    }
    mt->SetNumberOfThreads(kRequiredThreads);

    rm->SetUserInitialization(new P001Detector(mat));
    rm->SetUserInitialization(new P001PhysicsList());
    rm->SetUserInitialization(new P001Actions());

    // Initializes geometry and EM physics tables only. No event loop is started.
    rm->Initialize();

    if(gP001GeneratedEvents.load()!=0) {
        delete rm;
        throw std::runtime_error("P001 generated events during deterministic initialization");
    }

    auto* gamma=G4Gamma::Gamma();
    G4EmCalculator calc;
    CsvTable table=readCsv(target);
    if(table.rows.empty()) {
        delete rm;
        throw std::runtime_error("P001 target table is empty");
    }

    fs::create_directories(outPath.parent_path());
    fs::path tmpPath=outPath;
    tmpPath += ".tmp";
    std::error_code ec;
    fs::remove(tmpPath,ec);

    std::ofstream out(tmpPath);
    if(!out) {
        delete rm;
        throw std::runtime_error("Cannot create temporary P001 output: "+tmpPath.string());
    }
    out<<std::setprecision(17);
    out<<"run_id,row_index,geometry_id,source_normalization_id,geant4_evaluation_energy_MeV,mc_mu_over_rho_cm2_g,generated_event_count\n";

    std::size_t rowsWritten=0;
    for(const auto& row:table.rows) {
        const int rowIndex=std::stoi(csvGet(table,row,"row_index"));
        const double energyMeV=std::stod(csvGet(table,row,"geant4_evaluation_energy_MeV"));
        if(!std::isfinite(energyMeV) || energyMeV<=0.0) {
            out.close(); fs::remove(tmpPath,ec); delete rm;
            throw std::runtime_error("Invalid P001 evaluation energy");
        }

        G4double muPerVolume=0.0;
        for(const std::string process: {"Rayl","phot","compt","conv"}) {
            const G4double xs=calc.ComputeCrossSectionPerVolume(
                energyMeV*MeV,gamma,process,mat
            );
            if(!std::isfinite(xs) || xs<0.0) {
                out.close(); fs::remove(tmpPath,ec); delete rm;
                throw std::runtime_error(
                    "Invalid P001 Geant4 cross section for process "+process
                );
            }
            muPerVolume += xs;
        }

        const double muCmInverse=muPerVolume/(1.0/cm);
        const double densityGPerCm3=mat->GetDensity()/(g/cm3);
        const double muOverRho=muCmInverse/densityGPerCm3;

        if(!std::isfinite(muOverRho) || muOverRho<=0.0) {
            out.close(); fs::remove(tmpPath,ec); delete rm;
            throw std::runtime_error("Invalid P001 mass attenuation coefficient");
        }

        out<<"P001_G4_EM_COEFFICIENTS,"
           <<rowIndex
           <<",P001_INFINITE_MEDIUM_COEFFICIENT,NOT_APPLICABLE,"
           <<energyMeV<<","
           <<muOverRho<<","
           <<gP001GeneratedEvents.load()
           <<"\n";
        ++rowsWritten;
    }

    out.flush();
    if(!out) {
        out.close(); fs::remove(tmpPath,ec); delete rm;
        throw std::runtime_error("Failed while writing P001 output");
    }
    out.close();

    if(rowsWritten!=table.rows.size()) {
        fs::remove(tmpPath,ec); delete rm;
        throw std::runtime_error("P001 output row count differs from target row count");
    }
    if(gP001GeneratedEvents.load()!=0) {
        fs::remove(tmpPath,ec); delete rm;
        throw std::runtime_error("P001 transported events; deterministic result rejected");
    }

    fs::remove(outPath,ec);
    ec.clear();
    fs::rename(tmpPath,outPath,ec);
    if(ec) {
        fs::remove(tmpPath,ec);
        delete rm;
        throw std::runtime_error("Could not atomically install P001 output");
    }

    delete rm;
    return 0;
#endif
}

int main(int argc,char** argv) {
    try {
        if(argc<2){std::cerr<<"usage: phase1_geant4_runner info | p001 <target.csv> <out.csv> | neutron <config.ini> | photon_pdd <config.ini>\n";return 2;}
        std::string mode=argv[1];
        if(mode=="info") {
            std::cout<<"IMPLEMENTATION_LANGUAGE=C++17\nREQUIRED_THREADS="<<kRequiredThreads<<"\n";
#ifdef G4MULTITHREADED
            std::cout<<"MULTITHREADED=1\n";
#else
            std::cout<<"MULTITHREADED=0\n";
#endif
            std::cout<<"GEANT4_VERSION="<<G4Version<<"\n"; return 0;
        }
        if(mode=="p001" && argc==4) return runP001(argv[2],argv[3]);
        if(mode=="neutron" && argc==3) return runNeutron(argv[2]);
        if(mode=="photon_pdd" && argc==3) return runPhotonPdd(argv[2]);
        throw std::runtime_error("Invalid arguments");
    } catch(const std::exception& e){std::cerr<<"FATAL: "<<e.what()<<"\n";return 1;}
}
